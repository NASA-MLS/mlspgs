! Copyright (c 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

MODULE PVM ! Interface to the f77 pvm library.

  ! This module is a low level interface between the f77 pvm library and the
  ! MLS code (mainly aimed at level 2

  USE MLSCommon, ONLY: r8

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130), PRIVATE :: Id = & 
       "$Id$"
  CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! First we define constants as enumerated types effectively.  This is
  ! basically taken straight from the fpvm3.h file that comes with pvm, and as
  ! such will need to be updated if a newer version of pvm is installed.

  ! --------------------
  ! spawn 'flag' options
  ! --------------------
  INTEGER, PARAMETER :: PVMTASKDEFAULT    =  0
  INTEGER, PARAMETER :: PVMTASKHOST       =  1
  INTEGER, PARAMETER :: PVMTASKARCH       =  2
  INTEGER, PARAMETER :: PVMTASKDEBUG      =  4
  INTEGER, PARAMETER :: PVMTASKTRACE      =  8
  INTEGER, PARAMETER :: PVMMPPFRONT       = 16
  INTEGER, PARAMETER :: PVMHOSTCOMPL      = 32
  INTEGER, PARAMETER :: PVMNOSPAWNPARENT  = 64

  ! --------------------------------
  ! old option names still supported
  ! --------------------------------
  INTEGER, PARAMETER :: PVMHOST  =  1
  INTEGER, PARAMETER :: PVMARCH  =  2
  INTEGER, PARAMETER :: PVMDEBUG =  4
  INTEGER, PARAMETER :: PVMTRACE =  8

  ! -------------------------
  ! buffer 'encoding' options
  ! -------------------------
  INTEGER, PARAMETER ::  PVMDATADEFAULT = 0
  INTEGER, PARAMETER ::  PVMDATARAW     = 1
  INTEGER, PARAMETER ::  PVMDATAINPLACE = 2
  INTEGER, PARAMETER ::  PVMDATATRACE   = 4

  ! --------------------------------
  ! old option names still supported
  ! --------------------------------
  INTEGER, PARAMETER :: PVMDEFAULT = 0
  INTEGER, PARAMETER :: PVMRAW     = 1
  INTEGER, PARAMETER :: PVMINPLACE = 2

  ! ----------------------
  ! notify 'about' options
  ! ----------------------
  INTEGER, PARAMETER :: PVMTASKEXIT     = 1 
  INTEGER, PARAMETER :: PVMHOSTDELETE   = 2 
  INTEGER, PARAMETER :: PVMHOSTADD      = 3 
  INTEGER, PARAMETER :: PVMROUTEADD     = 4 
  INTEGER, PARAMETER :: PVMROUTEDELETE  = 5 
  INTEGER, PARAMETER :: PVMNOTIFYCANCEL = 256 

  ! --------------------------------
  ! packing/unpacking 'what' options
  ! --------------------------------
  INTEGER, PARAMETER :: PVMSTRING   = 0
  INTEGER, PARAMETER :: PVMBYTE1    = 1
  INTEGER, PARAMETER :: PVMINTEGER2 = 2
  INTEGER, PARAMETER :: PVMINTEGER4 = 3
  INTEGER, PARAMETER :: PVMREAL4    = 4
  INTEGER, PARAMETER :: PVMCOMPLEX8 = 5
  INTEGER, PARAMETER :: PVMREAL8    = 6
  INTEGER, PARAMETER :: PVMCOMPLEX16= 7

  ! --------------------------------
  ! setopt/getopt options for 'what'
  ! --------------------------------
  INTEGER, PARAMETER :: VMROUTE          = 1
  INTEGER, PARAMETER :: PVMDEBUGMASK     = 2
  INTEGER, PARAMETER :: PVMAUTOERR       = 3
  INTEGER, PARAMETER :: PVMOUTPUTTID     = 4
  INTEGER, PARAMETER :: PVMOUTPUTCODE    = 5
  INTEGER, PARAMETER :: PVMTRACETID      = 6
  INTEGER, PARAMETER :: PVMTRACECODE     = 7
  INTEGER, PARAMETER :: PVMTRACEBUFFER   = 8
  INTEGER, PARAMETER :: PVMTRACEOPTIONS  = 9
  INTEGER, PARAMETER :: PVMFRAGSIZE      = 10
  INTEGER, PARAMETER :: PVMRESVTIDS      = 11
  INTEGER, PARAMETER :: PVMSOUTPUTTID    = 12
  INTEGER, PARAMETER :: PVMSOUTPUTCODE   = 13
  INTEGER, PARAMETER :: PVMSTRACETID     = 14
  INTEGER, PARAMETER :: PVMSTRACECODE    = 15
  INTEGER, PARAMETER :: PVMSTRACEBUFFER  = 16
  INTEGER, PARAMETER :: PVMSTRACEOPTIONS = 17
  INTEGER, PARAMETER :: PVMSHOWTIDS      = 18
  INTEGER, PARAMETER :: PVMPOLLTYPE      = 19
  INTEGER, PARAMETER :: PVMPOLLTIME      = 20
  INTEGER, PARAMETER :: PVMOUTPUTCTX     = 21
  INTEGER, PARAMETER :: PVMTRACECTX      = 22
  INTEGER, PARAMETER :: PVMSOUTPUTCTX    = 23
  INTEGER, PARAMETER :: PVMSTRACECTX     = 24
  INTEGER, PARAMETER :: PVMNORESET       = 25

  ! --------------------------------------------
  ! tracing option values for setopt function
  ! --------------------------------------------
  INTEGER, PARAMETER :: PVMTRACEFULL     = 1
  INTEGER, PARAMETER :: PVMTRACETIME     = 2
  INTEGER, PARAMETER :: PVMTRACECOUNT    = 3

  ! --------------------------------------------
  ! poll type options for 'how' in setopt function
  ! --------------------------------------------
  INTEGER, PARAMETER :: PVMPOLLCONSTANT = 1
  INTEGER, PARAMETER :: PVMPOLLSLEEP    = 2

  ! --------------------------------------------
  ! for message mailbox operations
  ! --------------------------------------------
  INTEGER, PARAMETER :: PVMMBOXDEFAULT        =  0
  INTEGER, PARAMETER :: PVMMBOXPERSISTENT     =  1
  INTEGER, PARAMETER :: PVMMBOXMULTIINSTANCE  =  2
  INTEGER, PARAMETER :: PVMMBOXOVERWRITABLE   =  4
  INTEGER, PARAMETER :: PVMMBOXFIRSTAVAIL     =  8
  INTEGER, PARAMETER :: PVMMBOXREADANDDELETE  = 16
  INTEGER, PARAMETER :: PVMMBOXWAITFORINFO    = 32

  ! --------------------------------------------
  ! routing options for 'how' in setopt function
  ! --------------------------------------------
  INTEGER, PARAMETER :: PVMDONTROUTE  = 1
  INTEGER, PARAMETER :: PVMALLOWDIRECT= 2
  INTEGER, PARAMETER :: PVMROUTEDIRECT= 3

  ! --------------------------
  ! error 'info' return values
  ! --------------------------
  INTEGER, PARAMETER :: PvmOk           =   0
  INTEGER, PARAMETER :: PvmBadParam     =  -2
  INTEGER, PARAMETER :: PvmMismatch     =  -3
  INTEGER, PARAMETER :: PvmOverflow     =  -4
  INTEGER, PARAMETER :: PvmNoData       =  -5
  INTEGER, PARAMETER :: PvmNoHost       =  -6
  INTEGER, PARAMETER :: PvmNoFile       =  -7
  INTEGER, PARAMETER :: PvmDenied       =  -8
  INTEGER, PARAMETER :: PvmNoMem        = -10
  INTEGER, PARAMETER :: PvmBadMsg       = -12
  INTEGER, PARAMETER :: PvmSysErr       = -14
  INTEGER, PARAMETER :: PvmNoBuf        = -15
  INTEGER, PARAMETER :: PvmNoSuchBuf    = -16
  INTEGER, PARAMETER :: PvmNullGroup    = -17
  INTEGER, PARAMETER :: PvmDupGroup     = -18
  INTEGER, PARAMETER :: PvmNoGroup      = -19
  INTEGER, PARAMETER :: PvmNotInGroup   = -20
  INTEGER, PARAMETER :: PvmNoInst       = -21
  INTEGER, PARAMETER :: PvmHostFail     = -22
  INTEGER, PARAMETER :: PvmNoParent     = -23
  INTEGER, PARAMETER :: PvmNotImpl      = -24
  INTEGER, PARAMETER :: PvmDSysErr      = -25
  INTEGER, PARAMETER :: PvmBadVersion   = -26
  INTEGER, PARAMETER :: PvmOutOfRes     = -27
  INTEGER, PARAMETER :: PvmDupHost      = -28
  INTEGER, PARAMETER :: PvmCantStart    = -29
  INTEGER, PARAMETER :: PvmAlready      = -30
  INTEGER, PARAMETER :: PvmNoTask       = -31
  INTEGER, PARAMETER :: PvmNotFound     = -32
  INTEGER, PARAMETER :: PvmExists       = -33
  INTEGER, PARAMETER :: PvmHostrNMstr   = -34
  INTEGER, PARAMETER :: PvmParentNotSet = -35

  ! --------------------------
  ! these are going away in the next version.
  ! use the replacements
  ! --------------------------
  INTEGER, PARAMETER :: PvmNoEntry    = -32
  INTEGER, PARAMETER :: PvmDupEntry   = -33

  INTERFACE

     SUBROUTINE pvmfmytid(tid)
       INTEGER, INTENT(OUT) :: tid
     END SUBROUTINE pvmfmytid

     SUBROUTINE pvmfinitsend(encoding, bufid)
       INTEGER, INTENT(IN) :: encoding
       INTEGER, INTENT(OUT) :: bufid
     END SUBROUTINE pvmfinitsend

     SUBROUTINE pvmfbcast(group, msgtag, info)
       CHARACTER (LEN=*), INTENT(IN) :: group
       INTEGER, INTENT(IN) :: msgtag
       INTEGER, INTENT(OUT) :: info
     END SUBROUTINE pvmfbcast

     SUBROUTINE pvmfsend(tid, msgtag, info)
       INTEGER, INTENT(IN) :: tid
       INTEGER, INTENT(IN) :: msgtag
       INTEGER, INTENT(OUT) :: info
     END SUBROUTINE pvmfsend

     SUBROUTINE pvmfnrecv(tid, msgtag, bufid)
       INTEGER, INTENT(IN) :: tid
       INTEGER, INTENT(IN) :: msgtag
       INTEGER, INTENT(OUT) :: bufid
     END SUBROUTINE pvmfnrecv
     
     INTEGER FUNCTION pvm_pkstr(line)
       CHARACTER (LEN=*) :: line
     END FUNCTION pvm_pkstr

     INTEGER FUNCTION pvm_pkint(values,num,stride)
       INTEGER :: values(*)
       INTEGER :: num
       INTEGER :: stride
     END FUNCTION pvm_pkint

     INTEGER FUNCTION pvm_pkdouble(values,num,stride)
       DOUBLE PRECISION :: values(*)
       INTEGER :: num
       INTEGER :: stride
     END FUNCTION pvm_pkdouble

     INTEGER FUNCTION pvm_upkstr(line)
       CHARACTER (LEN=*) :: line
     END FUNCTION pvm_upkstr

     INTEGER FUNCTION pvm_upkint(values,num,stride)
       INTEGER :: values(*)
       INTEGER :: num
       INTEGER :: stride
     END FUNCTION pvm_upkint

     INTEGER FUNCTION pvm_upkdouble(values,num,stride)
       DOUBLE PRECISION :: values(*)
       INTEGER :: num
       INTEGER :: stride
     END FUNCTION pvm_upkdouble

     SUBROUTINE pvmfgsize(group, gsize)
       CHARACTER (LEN=*), INTENT(IN) :: group
       INTEGER, INTENT(OUT) :: gsize
     END SUBROUTINE pvmfgsize

     SUBROUTINE pvmfjoingroup(group, inum)
       CHARACTER (LEN=*), INTENT(IN) :: group
       INTEGER, INTENT(OUT) :: inum
     END SUBROUTINE pvmfjoingroup

  END INTERFACE

  INTERFACE pvmf90pack
     MODULE PROCEDURE pvmf90packString, pvmf90packInteger, pvmf90packReal, &
          & pvmf90packIntarr1, pvmf90packIntarr2, pvmf90packIntarr3, &
          & pvmf90packRealarr1, pvmf90packRealarr2, pvmf90packRealarr3
  END INTERFACE

  INTERFACE pvmf90unpack
     MODULE PROCEDURE pvmf90unpackString, pvmf90unpackInteger, pvmf90unpackReal, &
          & pvmf90unpackIntarr1, pvmf90unpackIntarr2, pvmf90unpackIntarr3, &
          & pvmf90unpackRealarr1, pvmf90unpackRealarr2, pvmf90unpackRealarr3
  END INTERFACE

CONTAINS

  SUBROUTINE pvmf90packString(line,info)
    CHARACTER (LEN=*), INTENT(IN) :: line
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkstr(line)
  END SUBROUTINE pvmf90packString

  SUBROUTINE pvmf90packInteger(value,info)
    INTEGER, INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkint( (/value/) ,1,1)
  END SUBROUTINE pvmf90packInteger

  SUBROUTINE pvmf90packIntarr1(values,info)
    INTEGER, DIMENSION(:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    
    info=pvm_pkint(values,SIZE(values),1)
  END SUBROUTINE pvmf90packIntarr1

  SUBROUTINE pvmf90packIntarr2(values,info)
    INTEGER, DIMENSION(:,:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkint(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90packIntarr2

  SUBROUTINE pvmf90packIntarr3(values,info)
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkint(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90packIntarr3

  SUBROUTINE pvmf90packReal(value,info)
    REAL (r8), INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkdouble((/value/),1,1)
  END SUBROUTINE pvmf90packReal

  SUBROUTINE pvmf90packRealarr1(values,info)
    REAL (r8), DIMENSION(:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkdouble(values,SIZE(values),1)
  END SUBROUTINE pvmf90packRealarr1

  SUBROUTINE pvmf90packRealarr2(values,info)
    REAL (r8), DIMENSION(:,:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkdouble(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90packRealarr2

  SUBROUTINE pvmf90packRealarr3(values,info)
    REAL (r8), DIMENSION(:,:,:), INTENT(IN) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_pkdouble(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90packRealarr3

  ! ---------------------------------------------------------------------------

  SUBROUTINE pvmf90unpackString(line,info)
    CHARACTER (LEN=*), INTENT(OUT) :: line
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkstr(line)
  END SUBROUTINE pvmf90unpackString

  SUBROUTINE pvmf90unpackInteger(value,info)
    INTEGER, INTENT(OUT) :: value
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(1) :: tempValue
    info=pvm_upkint( tempValue ,1,1)
    value=tempValue(1)
  END SUBROUTINE pvmf90unpackInteger

  SUBROUTINE pvmf90unpackIntarr1(values,info)
    INTEGER, DIMENSION(:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    
    info=pvm_upkint(values,SIZE(values),1)
  END SUBROUTINE pvmf90unpackIntarr1

  SUBROUTINE pvmf90unpackIntarr2(values,info)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkint(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90unpackIntarr2

  SUBROUTINE pvmf90unpackIntarr3(values,info)
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkint(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90unpackIntarr3

  SUBROUTINE pvmf90unpackReal(value,info)
    REAL (r8), INTENT(OUT) :: value
    INTEGER, INTENT(OUT) :: info

    REAL (r8), DIMENSION(1) :: tempValue
    info=pvm_upkdouble(tempValue,1,1)
    value=tempValue(1)
  END SUBROUTINE pvmf90unpackReal

  SUBROUTINE pvmf90unpackRealarr1(values,info)
    REAL (r8), DIMENSION(:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkdouble(values,SIZE(values),1)
  END SUBROUTINE pvmf90unpackRealarr1

  SUBROUTINE pvmf90unpackRealarr2(values,info)
    REAL (r8), DIMENSION(:,:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkdouble(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90unpackRealarr2

  SUBROUTINE pvmf90unpackRealarr3(values,info)
    REAL (r8), DIMENSION(:,:,:), INTENT(OUT) :: values
    INTEGER, INTENT(OUT) :: info
    info=pvm_upkdouble(RESHAPE(values,(/SIZE(values)/)),SIZE(values),1)
  END SUBROUTINE pvmf90unpackRealarr3

END MODULE PVM
