! Copyright (c 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PVM ! Interface to the f77 pvm library.

  ! This module is a low level interface between the f77 pvm library and the
  ! MLS code (mainly aimed at level 2

  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

  implicit none

  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130), private :: Id = & 
       "$Id$"
  character(LEN=*), private, parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! First we define constants as enumerated types effectively.  This is
  ! basically taken straight from the fpvm3.h file that comes with pvm, and as
  ! such will need to be updated if a newer version of pvm is installed.

  ! --------------------
  ! spawn 'flag' options
  ! --------------------
  integer, parameter :: PVMTASKDEFAULT    =  0
  integer, parameter :: PVMTASKHOST       =  1
  integer, parameter :: PVMTASKARCH       =  2
  integer, parameter :: PVMTASKDEBUG      =  4
  integer, parameter :: PVMTASKTRACE      =  8
  integer, parameter :: PVMMPPFRONT       = 16
  integer, parameter :: PVMHOSTCOMPL      = 32
  integer, parameter :: PVMNOSPAWNPARENT  = 64

  ! --------------------------------
  ! old option names still supported
  ! --------------------------------
  integer, parameter :: PVMHOST  =  1
  integer, parameter :: PVMARCH  =  2
  integer, parameter :: PVMDEBUG =  4
  integer, parameter :: PVMTRACE =  8

  ! -------------------------
  ! buffer 'encoding' options
  ! -------------------------
  integer, parameter ::  PVMDATADEFAULT = 0
  integer, parameter ::  PVMDATARAW     = 1
  integer, parameter ::  PVMDATAINPLACE = 2
  integer, parameter ::  PVMDATATRACE   = 4

  ! --------------------------------
  ! old option names still supported
  ! --------------------------------
  integer, parameter :: PVMDEFAULT = 0
  integer, parameter :: PVMRAW     = 1
  integer, parameter :: PVMINPLACE = 2

  ! ----------------------
  ! notify 'about' options
  ! ----------------------
  integer, parameter :: PVMTASKEXIT     = 1 
  integer, parameter :: PVMHOSTDELETE   = 2 
  integer, parameter :: PVMHOSTADD      = 3 
  integer, parameter :: PVMROUTEADD     = 4 
  integer, parameter :: PVMROUTEDELETE  = 5 
  integer, parameter :: PVMNOTIFYCANCEL = 256 

  ! --------------------------------
  ! packing/unpacking 'what' options
  ! --------------------------------
  integer, parameter :: PVMSTRING   = 0
  integer, parameter :: PVMBYTE1    = 1
  integer, parameter :: PVMINTEGER2 = 2
  integer, parameter :: PVMINTEGER4 = 3
  integer, parameter :: PVMREAL4    = 4
  integer, parameter :: PVMCOMPLEX8 = 5
  integer, parameter :: PVMREAL8    = 6
  integer, parameter :: PVMCOMPLEX16= 7

  ! --------------------------------
  ! setopt/getopt options for 'what'
  ! --------------------------------
  integer, parameter :: VMROUTE          = 1
  integer, parameter :: PVMDEBUGMASK     = 2
  integer, parameter :: PVMAUTOERR       = 3
  integer, parameter :: PVMOUTPUTTID     = 4
  integer, parameter :: PVMOUTPUTCODE    = 5
  integer, parameter :: PVMTRACETID      = 6
  integer, parameter :: PVMTRACECODE     = 7
  integer, parameter :: PVMTRACEBUFFER   = 8
  integer, parameter :: PVMTRACEOPTIONS  = 9
  integer, parameter :: PVMFRAGSIZE      = 10
  integer, parameter :: PVMRESVTIDS      = 11
  integer, parameter :: PVMSOUTPUTTID    = 12
  integer, parameter :: PVMSOUTPUTCODE   = 13
  integer, parameter :: PVMSTRACETID     = 14
  integer, parameter :: PVMSTRACECODE    = 15
  integer, parameter :: PVMSTRACEBUFFER  = 16
  integer, parameter :: PVMSTRACEOPTIONS = 17
  integer, parameter :: PVMSHOWTIDS      = 18
  integer, parameter :: PVMPOLLTYPE      = 19
  integer, parameter :: PVMPOLLTIME      = 20
  integer, parameter :: PVMOUTPUTCTX     = 21
  integer, parameter :: PVMTRACECTX      = 22
  integer, parameter :: PVMSOUTPUTCTX    = 23
  integer, parameter :: PVMSTRACECTX     = 24
  integer, parameter :: PVMNORESET       = 25

  ! --------------------------------------------
  ! tracing option values for setopt function
  ! --------------------------------------------
  integer, parameter :: PVMTRACEFULL     = 1
  integer, parameter :: PVMTRACETIME     = 2
  integer, parameter :: PVMTRACECOUNT    = 3

  ! --------------------------------------------
  ! poll type options for 'how' in setopt function
  ! --------------------------------------------
  integer, parameter :: PVMPOLLCONSTANT = 1
  integer, parameter :: PVMPOLLSLEEP    = 2

  ! --------------------------------------------
  ! for message mailbox operations
  ! --------------------------------------------
  integer, parameter :: PVMMBOXDEFAULT        =  0
  integer, parameter :: PVMMBOXPERSISTENT     =  1
  integer, parameter :: PVMMBOXMULTIINSTANCE  =  2
  integer, parameter :: PVMMBOXOVERWRITABLE   =  4
  integer, parameter :: PVMMBOXFIRSTAVAIL     =  8
  integer, parameter :: PVMMBOXREADANDDELETE  = 16
  integer, parameter :: PVMMBOXWAITFORINFO    = 32

  ! --------------------------------------------
  ! routing options for 'how' in setopt function
  ! --------------------------------------------
  integer, parameter :: PVMDONTROUTE  = 1
  integer, parameter :: PVMALLOWDIRECT= 2
  integer, parameter :: PVMROUTEDIRECT= 3

  ! --------------------------
  ! error 'info' return values
  ! --------------------------
  integer, parameter :: PvmOk           =   0
  integer, parameter :: PvmBadParam     =  -2
  integer, parameter :: PvmMismatch     =  -3
  integer, parameter :: PvmOverflow     =  -4
  integer, parameter :: PvmNoData       =  -5
  integer, parameter :: PvmNoHost       =  -6
  integer, parameter :: PvmNoFile       =  -7
  integer, parameter :: PvmDenied       =  -8
  integer, parameter :: PvmNoMem        = -10
  integer, parameter :: PvmBadMsg       = -12
  integer, parameter :: PvmSysErr       = -14
  integer, parameter :: PvmNoBuf        = -15
  integer, parameter :: PvmNoSuchBuf    = -16
  integer, parameter :: PvmNullGroup    = -17
  integer, parameter :: PvmDupGroup     = -18
  integer, parameter :: PvmNoGroup      = -19
  integer, parameter :: PvmNotInGroup   = -20
  integer, parameter :: PvmNoInst       = -21
  integer, parameter :: PvmHostFail     = -22
  integer, parameter :: PvmNoParent     = -23
  integer, parameter :: PvmNotImpl      = -24
  integer, parameter :: PvmDSysErr      = -25
  integer, parameter :: PvmBadVersion   = -26
  integer, parameter :: PvmOutOfRes     = -27
  integer, parameter :: PvmDupHost      = -28
  integer, parameter :: PvmCantStart    = -29
  integer, parameter :: PvmAlready      = -30
  integer, parameter :: PvmNoTask       = -31
  integer, parameter :: PvmNotFound     = -32
  integer, parameter :: PvmExists       = -33
  integer, parameter :: PvmHostrNMstr   = -34
  integer, parameter :: PvmParentNotSet = -35

  ! --------------------------
  ! these are going away in the next version.
  ! use the replacements
  ! --------------------------
  integer, parameter :: PvmNoEntry    = -32
  integer, parameter :: PvmDupEntry   = -33

  interface

     subroutine pvmfspawn ( task, flag, whr, ntask, tids, numt )
       character(len=*), intent(in) :: task
       integer, intent(in) :: flag
       character(len=*), intent(in) :: whr
       integer, intent(in) :: ntask
       integer, dimension(ntask), intent(out) :: tids
       integer, intent(out) :: numt
     end subroutine pvmfspawn

     subroutine pvmfmcast( ntask, tids, msgtag, info )
       integer,intent(in)::ntask,msgtag
       integer,intent(in),dimension(ntask)::tids
       integer,intent(out)::info
     end subroutine pvmfmcast

     subroutine pvmfcatchout ( onoff, info )
       integer, intent(in) :: ONOFF
       integer, intent(out) :: inFO
     end subroutine pvmfcatchout

     subroutine pvmfmytid(tid)
       integer, intent(out) :: tid
     end subroutine pvmfmytid

     subroutine pvmfinitsend(encoding, bufid)
       integer, intent(in) :: encoding
       integer, intent(out) :: bufid
     end subroutine pvmfinitsend

     subroutine pvmfbcast(group, msgtag, info)
       character (len=*), intent(in) :: group
       integer, intent(in) :: msgtag
       integer, intent(out) :: info
     end subroutine pvmfbcast

     subroutine pvmfbufinfo(bufid, bytes, msgtag, tid, info)
       integer, intent(in) :: bufID
       integer, intent(out) :: bytes
       integer, intent(out) :: msgtag
       integer, intent(out) :: tid
       integer, intent(out) :: info
     end subroutine pvmfbufinfo

     subroutine pvmfsend(tid, msgtag, info)
       integer, intent(in) :: tid
       integer, intent(in) :: msgtag
       integer, intent(out) :: info
     end subroutine pvmfsend

     subroutine pvmfnrecv(tid, msgtag, bufid)
       integer, intent(in) :: tid
       integer, intent(in) :: msgtag
       integer, intent(out) :: bufid
     end subroutine pvmfnrecv
     
     subroutine pvmfrecv(tid, msgtag, bufid)
       integer, intent(in) :: tid
       integer, intent(in) :: msgtag
       integer, intent(out) :: bufid
     end subroutine pvmfrecv
     
     integer function pvm_pkstr(line)
       character(len=*), intent(in) :: line
     end function pvm_pkstr

     integer function pvm_pkbyte(values,num,stride)
       character(len=1), intent(in) :: values(*)
       integer, intent(in) :: num
       integer, intent(in) :: stride
     end function pvm_pkbyte

     integer function pvm_pkint(values,num,stride)
       integer, intent(in) :: values(*)
       integer, intent(in) :: num
       integer, intent(in) :: stride
     end function pvm_pkint

     integer function pvm_pkdouble(values,num,stride)
       double precision, intent(in) :: values(*)
       integer, intent(in) :: num
       integer, intent(in) :: stride
     end function pvm_pkdouble

     integer function pvm_upkstr(line)
       character(len=*), intent(out) :: line
     end function pvm_upkstr

     integer function pvm_upkbyte(values,num,stride)
       character(len=1), intent(out) :: values(*)
       integer, intent(in) :: num
       integer, intent(in) :: stride
     end function pvm_upkbyte

     integer function pvm_upkint(values,num,stride)
       integer, intent(out) :: values(*)
       integer, intent(in) :: num
       integer, intent(in) :: stride
     end function pvm_upkint

     integer function pvm_upkdouble(values,num,stride)
       double precision, intent(out) :: values(*)
       integer, intent(in) :: num
       integer, intent(in)  :: stride
     end function pvm_upkdouble

     subroutine pvmfgsize(group, gsize)
       character (len=*), intent(in) :: group
       integer, intent(out) :: gsize
     end subroutine pvmfgsize

     subroutine pvmfjoingroup(group, inum)
       character (len=*), intent(in) :: group
       integer, intent(out) :: inum
     end subroutine pvmfjoingroup

     subroutine pvmfnotify(what, msgtag, cnt, tids, info)
       integer, intent(in) :: what
       integer, intent(in) :: msgtag
       integer, intent(in) :: cnt
       integer, intent(in),dimension(cnt) :: tids
       integer, intent(out) :: info
     end subroutine pvmfnotify

     ! These ones are to get around the irritating inability of pvmfspawn
     ! to pass arguments

     subroutine ClearPVMArgs()
     end subroutine ClearPVMArgs

     subroutine FreePVMArgs()
     end subroutine FreePVMArgs
     
     subroutine NextPVMArg(arg)
       character (len=*), intent(in) :: arg
     end subroutine NextPVMArg
     
     integer function MyPVMSpawn ( task, flag, whr, ntask, tids )
       character (len=*), intent(in) :: task
       integer, intent(in) :: flag
       character (len=*), intent(in) :: whr
       integer, intent(in) :: ntask
       integer, dimension(ntask), intent(out) :: tids
     end function MyPVMSpawn

  end interface

  interface pvmf90pack
     module procedure pvmf90packString, pvmf90packInteger, &
       & pvmf90packReal, pvmf90packChararr1, pvmf90packChararr2, &
       & pvmf90packIntarr1, pvmf90packIntarr2, pvmf90packIntarr3, &
       & pvmf90packRealarr1, pvmf90packRealarr2, pvmf90packRealarr3
  end interface

  interface pvmf90unpack
     module procedure pvmf90unpackString, pvmf90unpackInteger, &
       & pvmf90unpackReal, pvmf90unpackChararr1, pvmf90unpackChararr2, &
       & pvmf90unpackIntarr1, pvmf90unpackIntarr2, pvmf90unpackIntarr3, &
       & pvmf90unpackRealarr1, pvmf90unpackRealarr2, pvmf90unpackRealarr3
  end interface

contains

  subroutine pvmf90packString(line,info)
    character (LEN=*), intent(in) :: line
    integer, intent(out) :: info
    info=pvm_pkstr(line)
  end subroutine pvmf90packString

  subroutine pvmf90packInteger(value,info)
    integer, intent(in) :: value
    integer, intent(out) :: info
    info=pvm_pkint( (/value/) ,1,1)
  end subroutine pvmf90packInteger

  subroutine pvmf90packChararr1(values,info)
    character(len=1), dimension(:), intent(in) :: values
    integer, intent(out) :: info
    
    info=pvm_pkbyte(values,size(values),1)
  end subroutine pvmf90packChararr1

  subroutine pvmf90packChararr2(values,info)
    character(len=1), dimension(:,:), intent(in) :: values
    integer, intent(out) :: info
    
    info=pvm_pkbyte(values,size(values),1)
  end subroutine pvmf90packChararr2

  subroutine pvmf90packIntarr1(values,info)
    integer, dimension(:), intent(in) :: values
    integer, intent(out) :: info
    
    info=pvm_pkint(values,size(values),1)
  end subroutine pvmf90packIntarr1

  subroutine pvmf90packIntarr2(values,info)
    integer, dimension(:,:), intent(in) :: values
    integer, intent(out) :: info
    info=pvm_pkint(reshape(values,(/size(values)/)),size(values),1)
  end subroutine pvmf90packIntarr2

  subroutine pvmf90packIntarr3(values,info)
    integer, dimension(:,:,:), intent(in) :: values
    integer, intent(out) :: info
    info=pvm_pkint(reshape(values,(/size(values)/)),size(values),1)
  end subroutine pvmf90packIntarr3

  subroutine pvmf90packReal(value,info)
    real (r8), intent(in) :: value
    integer, intent(out) :: info
    info=pvm_pkdouble((/value/),1,1)
  end subroutine pvmf90packReal

  subroutine pvmf90packRealarr1(values,info)
    real (r8), dimension(:), intent(in) :: values
    integer, intent(out) :: info
    info=pvm_pkdouble(values,size(values),1)
  end subroutine pvmf90packRealarr1

  subroutine pvmf90packRealarr2(values,info)
    real (r8), dimension(:,:), intent(in) :: values
    integer, intent(out) :: info
    info=pvm_pkdouble(reshape(values,(/size(values)/)),size(values),1)
  end subroutine pvmf90packRealarr2

  subroutine pvmf90packRealarr3(values,info)
    real (r8), dimension(:,:,:), intent(in) :: values
    integer, intent(out) :: info
    info=pvm_pkdouble(reshape(values,(/size(values)/)),size(values),1)
  end subroutine pvmf90packRealarr3

  ! ---------------------------------------------------------------------------

  subroutine pvmf90unpackString(line,info)
    character (LEN=*), intent(out) :: line
    integer, intent(out) :: info
    info=pvm_upkstr(line)
  end subroutine pvmf90unpackString

  subroutine pvmf90unpackInteger(value,info)
    integer, intent(out) :: value
    integer, intent(out) :: info

    integer, dimension(1) :: tempValue
    info=pvm_upkint( tempValue ,1,1)
    value=tempValue(1)
  end subroutine pvmf90unpackInteger

  subroutine pvmf90unpackChararr1(values,info)
    character(len=1), dimension(:), intent(out) :: values
    integer, intent(out) :: info
    
    info=pvm_upkbyte(values,size(values),1)
  end subroutine pvmf90unpackChararr1

  subroutine pvmf90unpackChararr2(values,info)
    character(len=1), dimension(:,:), intent(out) :: values
    integer, intent(out) :: info
    
    info=pvm_upkbyte(values,size(values),1)
  end subroutine pvmf90unpackChararr2

  subroutine pvmf90unpackIntarr1(values,info)
    integer, dimension(:), intent(out) :: values
    integer, intent(out) :: info
    
    info=pvm_upkint(values,size(values),1)
  end subroutine pvmf90unpackIntarr1

  subroutine pvmf90unpackIntarr2(values,info)
    integer, dimension(:,:), intent(out) :: values
    integer, intent(out) :: info

    integer, dimension(size(values)) :: tmpVal

    info=pvm_upkint( tmpVal, size(values), 1)
    values=reshape(tmpVal,shape(values))
  end subroutine pvmf90unpackIntarr2

  subroutine pvmf90unpackIntarr3(values,info)
    integer, dimension(:,:,:), intent(out) :: values
    integer, intent(out) :: info

    integer, dimension(size(values)) :: tmpVal

    info=pvm_upkint( tmpVal, size(values), 1)
    values=reshape(tmpVal,shape(values))
  end subroutine pvmf90unpackIntarr3

  subroutine pvmf90unpackReal(value,info)
    real (r8), intent(out) :: value
    integer, intent(out) :: info

    real (r8), dimension(1) :: tempValue
    info=pvm_upkdouble(tempValue,1,1)
    value=tempValue(1)
  end subroutine pvmf90unpackReal

  subroutine pvmf90unpackRealarr1(values,info)
    real (r8), dimension(:), intent(out) :: values
    integer, intent(out) :: info

    info=pvm_upkdouble( values, size(values), 1)
  end subroutine pvmf90unpackRealarr1

  subroutine pvmf90unpackRealarr2(values,info)
    real (r8), dimension(:,:), intent(out) :: values
    integer, intent(out) :: info

    double precision, dimension(size(values)) :: tmpVal

    info=pvm_upkdouble( tmpVal, size(values), 1)
    values=reshape(tmpVal,shape(values))
  end subroutine pvmf90unpackRealarr2

  subroutine pvmf90unpackRealarr3(values,info)
    real (r8), dimension(:,:,:), intent(out) :: values
    integer, intent(out) :: info

    double precision, dimension(size(values)) :: tmpVal

    info=pvm_upkdouble( tmpVal, size(values), 1)
    values=reshape(tmpVal,shape(values))
  end subroutine pvmf90unpackRealarr3

  ! --------------------------------------------  PVMERRORMESSAGE  -----
  subroutine PVMErrorMessage ( inFO, PLACE )
    ! This routine is called to log a PVM error
    integer, intent(in) :: INFO
    character (LEN=*) :: PLACE

    character (LEN=132) :: LINE

    write (line, * ) info
    call MLSMessage(MLSMSG_Error,ModuleName,'PVM error '//trim(place)//&
      ' Info='//trim(adjustl(line)))
  end subroutine PVMErrorMessage

end module PVM

! $Log$
! Revision 2.9  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.8  2002/02/01 23:49:45  livesey
! Added stuff for character arrays (i.e. not strings)
!
! Revision 2.7  2001/11/19 17:00:30  pumphrey
! Added interface for pvmfmcast
!
! Revision 2.6  2001/05/25 01:05:14  livesey
! Added pvmfnotify
!
! Revision 2.5  2001/05/24 19:37:47  livesey
! Embarassing bug fix in unpack arrays!
!
! Revision 2.4  2001/05/23 01:42:54  livesey
! Various new changes, wrappers etc.
!
! Revision 2.3  2001/03/15 05:21:52  livesey
! Added pvmfrecv
!
! Revision 2.2  2001/03/08 21:51:33  livesey
! Tidied stuff up
!
