! Copyright (c 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PVMIDL ! Communicate with and IDL (NJL's pvmlib) process using pvm.

  ! This module is an interface between the f90 MLSL2 code and an IDL process
  ! using the pvmlib IDL routines written by NJL.

  ! It allows for F90 and IDL to exchange integers and reals (r8) and arrays of
  ! the same upto 3D, and strings (not arrays of strings though as there are
  ! length issues.

  use MLSCommon, only : r8
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFMYTID, PVMF90PACK, PVMF90UNPACK

  implicit none
 
  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130), private :: Id = & 
       "$Id$"
  character(LEN=*), private, parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------


  interface PVMIDLpack
     module procedure PVMIDLpackstring, PVMIDLpackInteger, PVMIDLpackReal, &
          & PVMIDLpackIntarr1, PVMIDLpackIntarr2, PVMIDLpackIntarr3, &
          & PVMIDLpackRealarr1, PVMIDLpackRealarr2, PVMIDLpackRealarr3
  end interface

  interface PVMIDLunpack
     module procedure PVMIDLunpackstring, PVMIDLunpackInteger, PVMIDLunpackReal, &
          & PVMIDLunpackIntarr1, PVMIDLunpackIntarr2, PVMIDLunpackIntarr3, &
          & PVMIDLunpackRealarr1, PVMIDLunpackRealarr2, PVMIDLunpackRealarr3
  end interface

  interface PVMIDLSend
     module procedure PVMIDLSendString, PVMIDLSendInteger, PVMIDLSendReal,&
          & PVMIDLSendIntarr1, PVMIDLSendIntarr2, PVMIDLSendIntarr3, &
          & PVMIDLSendRealarr1, PVMIDLSendRealarr2, PVMIDLSendRealarr3
  end interface

  interface PVMIDLReceive
     module procedure PVMIDLReceiveString, PVMIDLReceiveInteger, PVMIDLReceiveReal,&
          & PVMIDLReceiveIntarr1, PVMIDLReceiveIntarr2, PVMIDLReceiveIntarr3, &
          & PVMIDLReceiveRealarr1, PVMIDLReceiveRealarr2, PVMIDLReceiveRealarr3
  end interface

  integer, parameter :: IDLMsgTag=100

contains

  subroutine PVMIDLpackString(line,info)
    character (LEN=*), intent(in) :: line
    integer, intent(out) :: info

    integer :: length

    ! First pack noDims and a 7 to indicate string
    call pvmf90pack( (/0,7/), info)

    ! Now pack the length of the string
    length=len_trim(line)
    if (info==0) call pvmf90pack( length, info)

    ! Now pack the string itself
    if ((info==0).and.(length/=0)) call pvmf90pack(trim(line),info)

  end subroutine PVMIDLpackString

  subroutine PVMIDLpackInteger(value,info)
    integer, intent(in) :: value
    integer, intent(out) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/0,3/), info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(value,info)

  end subroutine PVMIDLpackInteger

  subroutine PVMIDLpackReal(value,info)
    real (r8), intent(in) :: value
    integer, intent(out) :: info

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/0,5/), info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(value,info)
  end subroutine PVMIDLpackReal

  subroutine PVMIDLpackIntarr1(values,info)
    integer, intent(in), dimension(:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/1,3/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack( (/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackIntarr1

  subroutine PVMIDLpackIntarr2(values,info)
    integer, intent(in), dimension(:,:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/2,3/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack((/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackIntarr2

  subroutine PVMIDLpackIntarr3(values,info)
    integer, intent(in), dimension(:,:,:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/3,3/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack((/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackIntarr3

  subroutine PVMIDLpackRealarr1(values,info)
    real (r8), intent(in), dimension(:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/1,5/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack((/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackRealarr1

  subroutine PVMIDLpackRealarr2(values,info)
    real (r8), intent(in), dimension(:,:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/2,5/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack((/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackRealarr2

  subroutine PVMIDLpackRealarr3(values,info)
    real (r8), intent(in), dimension(:,:,:) :: values
    integer, intent(out) :: info

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/3,5/), info)

    ! Now output the dimensions themselves
    if (info==0) call pvmf90pack((/shape(values),size(values)/),info)

    ! Now pack the data itself
    if (info==0) call pvmf90pack(values,info)
  end subroutine PVMIDLpackRealarr3

  ! ------------------------------------------------------------------------------

  subroutine PVMIDLunpackString(line,info)
    character (LEN=*), intent(out) :: line
    integer, intent(out) :: info

    integer :: length
    integer, dimension(2) :: details

    ! First unpack noDims and a 7 to indicate string
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/=(/0,7/))) info= -200

       ! Now unpack the length of the string
       if (info==0) call pvmf90unpack( length, info)
       
       if ((info==0).and.(length > len(line))) info=-201
       
       ! Now unpack the string itself
       if ((info==0).and.(length/=0)) call pvmf90unpack(line,info)
       line=line(1:length)
    endif

  end subroutine PVMIDLunpackString

  subroutine PVMIDLunpackInteger(value,info)
    integer, intent(out) :: value
    integer, intent(out) :: info

    integer, dimension(2) :: details

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/=(/0,3/))) info= -200

       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(value,info)
    end if
  end subroutine PVMIDLunpackInteger

  subroutine PVMIDLunpackReal(value,info)
    real (r8), intent(out) :: value
    integer, intent(out) :: info

    integer, dimension(2) :: details

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/0,5/)) ) info= -200

       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(value,info)
    end if
  end subroutine PVMIDLunpackReal
     
  subroutine PVMIDLunpackIntarr1(values,info)
    integer, intent(out), dimension(:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(2) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/1,3/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:1)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackIntarr1

  subroutine PVMIDLunpackIntarr2(values,info)
    integer, intent(out), dimension(:,:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(3) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/2,3/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:2)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackIntarr2

  subroutine PVMIDLunpackIntarr3(values,info)
    integer, intent(out), dimension(:,:,:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(3) :: sentShape

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/3,3/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:3)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackIntarr3
     
  subroutine PVMIDLunpackRealarr1(values,info)
    real (r8), intent(out), dimension(:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(2) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/1,5/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:1)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackRealarr1

  subroutine PVMIDLunpackRealarr2(values,info)
    real (r8), intent(out), dimension(:,:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(3) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/2,5/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:2)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackRealarr2

  subroutine PVMIDLunpackRealarr3(values,info)
    real (r8), intent(out), dimension(:,:,:) :: values
    integer, intent(out) :: info

    integer, dimension(2) :: details
    integer, dimension(3) :: sentShape

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, info)

    if (info==0) then 
       if (any(details/= (/3,5/)) ) info= -200

       ! Now output the dimensions themselves
       if (info==0) call pvmf90unpack( sentShape,info)
       if (any(sentShape(1:3)/=shape(values))) info= -201
       
       ! Now unpack the data itself
       if (info==0) call pvmf90unpack(values,info)
    end if
  end subroutine PVMIDLunpackRealarr3

  ! ----------------------------------------------------------------------

  ! Now some simple all in one routines to do the initsend, the pack and the
  ! send all at once.  This should make life easier

  subroutine PVMIDLSendString(value,tid,info,msgTag)
    character (LEN=*), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, MYMSGTAG

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid,myMsgTag,info)
  end subroutine PVMIDLSendString

  subroutine PVMIDLSendInteger(value,tid,info,msgTag)
    integer, intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, MYMSGTAG

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value, info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendInteger

  subroutine PVMIDLSendReal(value,tid, info, msgTag)
    real(r8), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendReal


  subroutine PVMIDLSendIntArr1(value,tid, info, msgTag)
    integer, dimension(:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendIntArr1

  subroutine PVMIDLSendIntArr2(value,tid, info, msgTag)
    integer, dimension(:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendIntArr2

  subroutine PVMIDLSendIntArr3(value,tid, info, msgTag)
    integer, dimension(:,:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendIntArr3

  subroutine PVMIDLSendRealArr1(value,tid, info, msgTag)
    real(r8), dimension(:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag, info)
  end subroutine PVMIDLSendRealArr1

  subroutine PVMIDLSendRealArr2(value,tid, info, msgTag)
    real(r8), dimension(:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag,info)
  end subroutine PVMIDLSendRealArr2

  subroutine PVMIDLSendRealArr3(value,tid, info, msgTag)
    real(r8), dimension(:,:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    integer, intent(in), optional :: msgTag

    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,info)
    if (info==0) call PVMFSend(tid, myMsgTag, info)
  end subroutine PVMIDLSendRealArr3

  ! ----------------------------------------------------------------------

  ! Now the same for receive

  subroutine PVMIDLReceiveString(value,tid,info, noBlock, msgTag)
    character (LEN=*), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveString

  subroutine PVMIDLReceiveInteger(value,tid,info, noBlock, msgTag)
    integer, intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveInteger

  subroutine PVMIDLReceiveReal(value,tid,info, noBlock, msgTag)
    real(r8), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveReal


  subroutine PVMIDLReceiveIntArr1(value,tid,info, noBlock, msgTag)
    integer, dimension(:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveIntArr1

  subroutine PVMIDLReceiveIntArr2(value,tid,info, noBlock, msgTag)
    integer, dimension(:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveIntArr2

  subroutine PVMIDLReceiveIntArr3(value,tid,info, noBlock, msgTag)
    integer, dimension(:,:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveIntArr3

  subroutine PVMIDLReceiveRealArr1(value,tid,info, noBlock, msgTag)
    real(r8), dimension(:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveRealArr1

  subroutine PVMIDLReceiveRealArr2(value,tid,info, noBlock, msgTag)
    real(r8), dimension(:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveRealArr2

  subroutine PVMIDLReceiveRealArr3(value,tid,info, noBlock, msgTag)
    real(r8), dimension(:,:,:), intent(out) :: value
    integer, intent(in):: tid
    integer, intent(out) :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (present(noBlock)) useNoBlock=noBlock
    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info)
  end subroutine PVMIDLReceiveRealArr3

end module PVMIDL



