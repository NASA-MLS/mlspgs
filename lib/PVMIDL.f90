! Copyright (c 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

MODULE PVMIDL ! Communicate with and IDL (NJL's pvmlib) process using pvm.

  ! This module is an interface between the f90 MLSL2 code and an IDL process
  ! using the pvmlib IDL routines written by NJL.

  ! It allows for F90 and IDL to exchange integers and reals (r8) and arrays of
  ! the same upto 3D, and strings (not arrays of strings though as there are
  ! length issues.

  USE MLSCommon, ONLY : r8
  USE PVM

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130), PRIVATE :: Id = & 
       "$Id$"
  CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTERFACE PVMIDLpack
     MODULE PROCEDURE PVMIDLpackstring, PVMIDLpackInteger, PVMIDLpackReal, &
          & PVMIDLpackIntarr1, PVMIDLpackIntarr2, PVMIDLpackIntarr3, &
          & PVMIDLpackRealarr1, PVMIDLpackRealarr2, PVMIDLpackRealarr3
  END INTERFACE

  INTERFACE PVMIDLunpack
     MODULE PROCEDURE PVMIDLunpackstring, PVMIDLunpackInteger, PVMIDLunpackReal, &
          & PVMIDLunpackIntarr1, PVMIDLunpackIntarr2, PVMIDLunpackIntarr3, &
          & PVMIDLunpackRealarr1, PVMIDLunpackRealarr2, PVMIDLunpackRealarr3
  END INTERFACE

  INTERFACE PVMIDLSend
     MODULE PROCEDURE PVMIDLSendString, PVMIDLSendInteger, PVMIDLSendReal,&
          & PVMIDLSendIntarr1, PVMIDLSendIntarr2, PVMIDLSendIntarr3, &
          & PVMIDLSendRealarr1, PVMIDLSendRealarr2, PVMIDLSendRealarr3
  END INTERFACE

  INTERFACE PVMIDLReceive
     MODULE PROCEDURE PVMIDLReceiveString, PVMIDLReceiveInteger, PVMIDLReceiveReal,&
          & PVMIDLReceiveIntarr1, PVMIDLReceiveIntarr2, PVMIDLReceiveIntarr3, &
          & PVMIDLReceiveRealarr1, PVMIDLReceiveRealarr2, PVMIDLReceiveRealarr3
  END INTERFACE

  INTEGER, PARAMETER :: IDLMsgTag=100

CONTAINS

  SUBROUTINE PVMIDLpackString(line,info)
    CHARACTER (LEN=*), INTENT(IN) :: line
    INTEGER, INTENT(OUT) :: info

    INTEGER :: length

    ! First pack noDims and a 7 to indicate string
    CALL pvmf90pack( (/0,7/), info)

    ! Now pack the length of the string
    length=LEN_TRIM(line)
    IF (info==0) CALL pvmf90pack( length, info)

    ! Now pack the string itself
    IF ((info==0).AND.(length/=0)) CALL pvmf90pack(TRIM(line),info)

  END SUBROUTINE PVMIDLpackString

  SUBROUTINE PVMIDLpackInteger(value,info)
    INTEGER, INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90pack( (/0,3/), info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(value,info)

  END SUBROUTINE PVMIDLpackInteger

  SUBROUTINE PVMIDLpackReal(value,info)
    REAL (r8), INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 5 to indicate double
    CALL pvmf90pack( (/0,5/), info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(value,info)
  END SUBROUTINE PVMIDLpackReal

  SUBROUTINE PVMIDLpackIntarr1(values,info)
    INTEGER, INTENT(IN), DIMENSION(:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90pack( (/1,3/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack( (/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackIntarr1

  SUBROUTINE PVMIDLpackIntarr2(values,info)
    INTEGER, INTENT(IN), DIMENSION(:,:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90pack( (/2,3/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack((/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackIntarr2

  SUBROUTINE PVMIDLpackIntarr3(values,info)
    INTEGER, INTENT(IN), DIMENSION(:,:,:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90pack( (/3,3/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack((/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackIntarr3

  SUBROUTINE PVMIDLpackRealarr1(values,info)
    REAL (r8), INTENT(IN), DIMENSION(:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 5 to indicate double
    CALL pvmf90pack( (/1,5/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack((/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackRealarr1

  SUBROUTINE PVMIDLpackRealarr2(values,info)
    REAL (r8), INTENT(IN), DIMENSION(:,:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 5 to indicate double
    CALL pvmf90pack( (/2,5/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack((/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackRealarr2

  SUBROUTINE PVMIDLpackRealarr3(values,info)
    REAL (r8), INTENT(IN), DIMENSION(:,:,:) :: values
    INTEGER, INTENT(OUT) :: info

    ! First pack noDims and a 5 to indicate double
    CALL pvmf90pack( (/3,5/), info)

    ! Now output the dimensions themselves
    IF (info==0) CALL pvmf90pack((/SHAPE(values),SIZE(values)/),info)

    ! Now pack the data itself
    IF (info==0) CALL pvmf90pack(values,info)
  END SUBROUTINE PVMIDLpackRealarr3

  ! ------------------------------------------------------------------------------

  SUBROUTINE PVMIDLunpackString(line,info)
    CHARACTER (LEN=*), INTENT(OUT) :: line
    INTEGER, INTENT(OUT) :: info

    INTEGER :: length
    INTEGER, DIMENSION(2) :: details

    ! First unpack noDims and a 7 to indicate string
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/=(/0,7/))) info= -200

       ! Now unpack the length of the string
       IF (info==0) CALL pvmf90unpack( length, info)
       
       IF ((info==0).AND.(length > LEN(line))) info=-201
       
       ! Now unpack the string itself
       IF ((info==0).AND.(length/=0)) CALL pvmf90unpack(line,info)
       line=line(1:length)
    ENDIF

  END SUBROUTINE PVMIDLunpackString

  SUBROUTINE PVMIDLunpackInteger(value,info)
    INTEGER, INTENT(OUT) :: value
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/=(/0,3/))) info= -200

       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(value,info)
    END IF
  END SUBROUTINE PVMIDLunpackInteger

  SUBROUTINE PVMIDLunpackReal(value,info)
    REAL (r8), INTENT(OUT) :: value
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details

    ! First unpack noDims and a 5 to indicate double
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/0,5/)) ) info= -200

       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(value,info)
    END IF
  END SUBROUTINE PVMIDLunpackReal
     
  SUBROUTINE PVMIDLunpackIntarr1(values,info)
    INTEGER, INTENT(OUT), DIMENSION(:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(2) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/1,3/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:1)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackIntarr1

  SUBROUTINE PVMIDLunpackIntarr2(values,info)
    INTEGER, INTENT(OUT), DIMENSION(:,:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(3) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/2,3/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:2)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackIntarr2

  SUBROUTINE PVMIDLunpackIntarr3(values,info)
    INTEGER, INTENT(OUT), DIMENSION(:,:,:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(3) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/3,3/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:3)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackIntarr3
     
  SUBROUTINE PVMIDLunpackRealarr1(values,info)
    REAL (r8), INTENT(OUT), DIMENSION(:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(2) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/1,5/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:1)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackRealarr1

  SUBROUTINE PVMIDLunpackRealarr2(values,info)
    REAL (r8), INTENT(OUT), DIMENSION(:,:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(3) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/2,5/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:2)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackRealarr2

  SUBROUTINE PVMIDLunpackRealarr3(values,info)
    REAL (r8), INTENT(OUT), DIMENSION(:,:,:) :: values
    INTEGER, INTENT(OUT) :: info

    INTEGER, DIMENSION(2) :: details
    INTEGER, DIMENSION(3) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    CALL pvmf90unpack( details, info)

    IF (info==0) THEN 
       IF (ANY(details/= (/3,5/)) ) info= -200

       ! Now output the dimensions themselves
       IF (info==0) CALL pvmf90unpack( sentShape,info)
       IF (ANY(sentShape(1:3)/=SHAPE(values))) info= -201
       
       ! Now unpack the data itself
       IF (info==0) CALL pvmf90unpack(values,info)
    END IF
  END SUBROUTINE PVMIDLunpackRealarr3

  ! ----------------------------------------------------------------------

  ! Now some simple all in one routines to do the initsend, the pack and the
  ! send all at once.  This should make life easier

  SUBROUTINE PVMIDLSendString(value,tid,info)
    CHARACTER (LEN=*), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendString

  SUBROUTINE PVMIDLSendInteger(value,tid,info)
    INTEGER, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendInteger

  SUBROUTINE PVMIDLSendReal(value,tid,info)
    REAL(r8), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendReal


  SUBROUTINE PVMIDLSendIntArr1(value,tid,info)
    INTEGER, DIMENSION(:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendIntArr1

  SUBROUTINE PVMIDLSendIntArr2(value,tid,info)
    INTEGER, DIMENSION(:,:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendIntArr2

  SUBROUTINE PVMIDLSendIntArr3(value,tid,info)
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendIntArr3

  SUBROUTINE PVMIDLSendRealArr1(value,tid,info)
    REAL(r8), DIMENSION(:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendRealArr1

  SUBROUTINE PVMIDLSendRealArr2(value,tid,info)
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendRealArr2

  SUBROUTINE PVMIDLSendRealArr3(value,tid,info)
    REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info

    INTEGER :: bufferID

    CALL PVMFInitSend(PvmDataDefault,bufferID)
    CALL PVMIDLPack(value,info)
    IF (info==0) CALL PVMFSend(tid,IDLMsgTag,info)
  END SUBROUTINE PVMIDLSendRealArr3

  ! ----------------------------------------------------------------------

  ! Now the same for receive

  SUBROUTINE PVMIDLReceiveString(value,tid,info,noBlock)
    CHARACTER (LEN=*), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveString

  SUBROUTINE PVMIDLReceiveInteger(value,tid,info,noBlock)
    INTEGER, INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveInteger

  SUBROUTINE PVMIDLReceiveReal(value,tid,info,noBlock)
    REAL(r8), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveReal


  SUBROUTINE PVMIDLReceiveIntArr1(value,tid,info,noBlock)
    INTEGER, DIMENSION(:), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveIntArr1

  SUBROUTINE PVMIDLReceiveIntArr2(value,tid,info,noBlock)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveIntArr2

  SUBROUTINE PVMIDLReceiveIntArr3(value,tid,info,noBlock)
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveIntArr3

  SUBROUTINE PVMIDLReceiveRealArr1(value,tid,info,noBlock)
    REAL(r8), DIMENSION(:), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveRealArr1

  SUBROUTINE PVMIDLReceiveRealArr2(value,tid,info,noBlock)
    REAL(r8), DIMENSION(:,:), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveRealArr2

  SUBROUTINE PVMIDLReceiveRealArr3(value,tid,info,noBlock)
    REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: value
    INTEGER, INTENT(IN):: tid
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: noBlock

    LOGICAL :: useNoBlock = .FALSE.
    INTEGER :: bufferID

    IF (PRESENT(noBlock)) useNoBlock=noBlock
    IF (useNoBlock) THEN
       CALL PVMFNrecv(tid,IDLMsgTag,bufferID)
    ELSE
       CALL PVMFrecv(tid,IDLMsgTag,bufferID)
    ENDIF       
    CALL PVMIDLUnpack(value,info)
  END SUBROUTINE PVMIDLReceiveRealArr3

END MODULE PVMIDL



