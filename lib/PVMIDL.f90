! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

module PVMIDL ! Communicate with and IDL (NJL's pvmlib) process using pvm.

  ! This module is an interface between the f90 MLSL2 code and an IDL process
  ! using the pvmlib IDL routines written by NJL.

  ! It allows for F90 and IDL to exchange integers, reals (r4,r8) and arrays of
  ! the same upto 3D, and strings (not arrays of strings though as there are
  ! length issues.

  use MLSCommon, only : r4, r8
  use MLSMessageModule, only : PVMErrorMessage
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFMYTID, &
    & PVMF90PACK, PVMF90UNPACK

  implicit none
  public
  private PVMIDLStat
 
  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)), private :: Id = idParm
  character(LEN=*), parameter, private :: ModuleName="$RCSfile$"
  private not_used_here
  !-----------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
! This module is an interface between the f90 MLSL2 code and an IDL process  
! using the pvmlib IDL routines written by NJL.                              
!     (subroutines and functions)
! PVMIDLpack                      
! PVMIDLunpack                    
! PVMIDLSend                      
! PVMIDLReceive                   
! === (end of toc) ===

! === (start of api) ===
! PVMIDLpack (arg, int info, [message] ) 
! PVMIDLunpack (arg, int info) 
! PVMIDLSend (arg, int tid, [log noBlock], [msgTag])
! PVMIDLReceive (arg, int tid, [log noBlock], [msgTag]) 
!     arg can be one of:
!    {char* line, int value, r4 value, r8 value, char line(:), char line(:,:),
!     int values(:), int values(:,:), int values(:,:,:),
!     r4 values(:), r4 values(:,:), r4 values(:,:,:),
!     r8 values(:), r8 values(:,:), r8 values(:,:,:)}
! === (end of api) ===

  interface PVMIDLpack
     module procedure PVMIDLpackstring, PVMIDLpackInteger, PVMIDLpackReal, &
          & PVMIDLpackLogical, PVMIDLpackChararr1, PVMIDLpackChararr2, &
          & PVMIDLpackIntarr1, PVMIDLpackIntarr2, PVMIDLpackIntarr3, &
          & PVMIDLpackRealarr1, PVMIDLpackRealarr2, PVMIDLpackRealarr3,&
          & PVMIDLpackLogArr1, PVMIDLpackSngl, &
          & PVMIDLpackSnglarr1, PVMIDLpackSnglarr2, PVMIDLpackSnglarr3
  end interface

  interface PVMIDLunpack
     module procedure PVMIDLunpackstring, PVMIDLunpackInteger, PVMIDLunpackReal, &
          & PVMIDLunpackLogical, PVMIDLunpackChararr1, PVMIDLunpackChararr2, &
          & PVMIDLunpackIntarr1, PVMIDLunpackIntarr2, PVMIDLunpackIntarr3, &
          & PVMIDLunpackRealarr1, PVMIDLunpackRealarr2, PVMIDLunpackRealarr3, &
          & PVMIDLunpackLogarr1, PVMIDLunpackSngl, &
          & PVMIDLunpackSnglarr1, PVMIDLunpackSnglarr2, PVMIDLunpackSnglarr3
  end interface

  interface PVMIDLSend
     module procedure PVMIDLSendString, PVMIDLSendInteger, PVMIDLSendReal, &
          & PVMIDLSendLogical, &
          & PVMIDLSendIntarr1, PVMIDLSendIntarr2, PVMIDLSendIntarr3, &
          & PVMIDLSendRealarr1, PVMIDLSendRealarr2, PVMIDLSendRealarr3, &
          & PVMIDLSendLogarr1
  end interface

  interface PVMIDLReceive
     module procedure PVMIDLReceiveString, PVMIDLReceiveInteger, PVMIDLReceiveReal, &
          & PVMIDLReceiveLogical, &
          & PVMIDLReceiveIntarr1, PVMIDLReceiveIntarr2, PVMIDLReceiveIntarr3, &
          & PVMIDLReceiveRealarr1, PVMIDLReceiveRealarr2, PVMIDLReceiveRealarr3, &
          & PVMIDLReceiveLogarr1
  end interface

  integer, parameter :: IDLMsgTag=100

contains

  subroutine PVMIDLpackString(line,info,msg)
    character (LEN=*), intent(in) :: line
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: length, stat

    ! First pack noDims and a 7 to indicate string
    call pvmf90pack( (/0,7/), stat)

    ! Now pack the length of the string
    length=len_trim(line)
    if (stat==0) call pvmf90pack( length, stat)

    ! Now pack the string itself
    if ((stat==0).and.(length/=0)) call pvmf90pack(trim(line),stat)

    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackString

  subroutine PVMIDLpackInteger(value,info,msg)
    integer, intent(in) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/0,3/), stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(value,info)
    call PVMIDLStat ( stat, info, msg )

  end subroutine PVMIDLpackInteger

  subroutine PVMIDLpackSngl(value,info,msg)
    real (r4), intent(in) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    real(r8) :: r8value
    integer :: stat

    r8value = value
    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/0,5/), stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(r8value,info)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackSngl

  subroutine PVMIDLpackReal(value,info,msg)
    real (r8), intent(in) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/0,5/), stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(value,info)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackReal

  subroutine PVMIDLpackLogical(value,info,msg)
    logical, intent(in) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    integer :: intValue
    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/0,3/), stat)

    ! Now pack the data itself
    intValue = 0
    if ( value ) intValue=1
    if (stat==0) call pvmf90pack(intValue,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackLogical

  subroutine PVMIDLpackChararr1 ( line,info,msg )
    character (LEN=1), dimension(:), intent(in) :: line
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: length, stat

    ! First pack noDims and a 1 to indicate byte
    call pvmf90pack( (/1,1/), stat)

    ! Now pack the length of the array
    length=size(line)
    if (stat==0) call pvmf90pack ( length, stat )

    ! Now pack the string itself
    if ((stat==0).and.(length/=0)) call pvmf90pack(line,stat)
    call PVMIDLStat ( stat, info, msg )

  end subroutine PVMIDLpackChararr1

  subroutine PVMIDLpackChararr2 ( line,info,msg )
    character (LEN=1), dimension(:,:), intent(in) :: line
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: length, stat

    ! First pack noDims and a 1 to indicate byte
    call pvmf90pack( (/2,1/), stat)

    length = size(line)
    ! Now pack the length of the array
    if (stat==0) call pvmf90pack ( (/shape(line),size(line)/), stat )

    ! Now pack the string itself
    if ((stat==0).and.(length/=0)) call pvmf90pack(line,stat)
    call PVMIDLStat ( stat, info, msg )

  end subroutine PVMIDLpackChararr2

  subroutine PVMIDLpackIntarr1 ( values,info,msg )
    integer, intent(in), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/1,3/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack( (/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,info)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackIntarr1

  subroutine PVMIDLpackIntarr2(values,info,msg)
    integer, intent(in), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/2,3/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackIntarr2

  subroutine PVMIDLpackIntarr3(values,info,msg)
    integer, intent(in), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/3,3/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackIntarr3

  subroutine PVMIDLpackRealarr1(values,info,msg)
    real (r8), intent(in), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/1,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackRealarr1

  subroutine PVMIDLpackRealarr2(values,info,msg)
    real (r8), intent(in), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/2,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackRealarr2

  subroutine PVMIDLpackRealarr3(values,info,msg)
    real (r8), intent(in), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/3,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(values,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackRealarr3

  subroutine PVMIDLpackSnglarr1(values,info,msg)
    real (r4), intent(in), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/1,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(dble(values),stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackSnglarr1

  subroutine PVMIDLpackSnglarr2(values,info,msg)
    real (r4), intent(in), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/2,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(dble(values),stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackSnglarr2

  subroutine PVMIDLpackSnglarr3(values,info,msg)
    real (r4), intent(in), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    ! First pack noDims and a 5 to indicate double
    call pvmf90pack( (/3,5/), stat)

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack((/shape(values),size(values)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(dble(values),stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackSnglarr3

  subroutine PVMIDLpackLogarr1(values,info,msg)
    logical, intent(in), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    integer :: stat

    integer, dimension(size(values)) :: valAsInt

    ! First pack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90pack( (/1,3/), stat)

    valAsInt = 0
    where (values)
      valAsInt = 1
    end where

    ! Now output the dimensions themselves
    if (stat==0) call pvmf90pack( (/shape(valAsInt),size(valAsInt)/),stat)

    ! Now pack the data itself
    if (stat==0) call pvmf90pack(valAsInt,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLpackLogarr1


  ! ------------------------------------------------------------------------------

  subroutine PVMIDLunpackString(line,info,msg)
    character (LEN=*), intent(out) :: line
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: length, stat
    integer, dimension(2) :: details

    ! First unpack noDims and a 7 to indicate string
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/=(/0,7/))) stat= -200

       ! Now unpack the length of the string
       if (stat==0) call pvmf90unpack( length, stat)
       
       if ((stat==0).and.(length > len(line))) stat=-201
       
       ! Now unpack the string itself
       if ((stat==0).and.(length/=0)) then
         call pvmf90unpack(line,stat)
         line = line(1:length)
       else
         line = ''
       end if
    end if
    call PVMIDLStat ( stat, info, msg )

  end subroutine PVMIDLunpackString

  subroutine PVMIDLunpackInteger(value,info,msg)
    integer, intent(out) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/=(/0,3/))) stat= -200

       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(value,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackInteger

  subroutine PVMIDLunpackReal(value,info,msg)
    real (r8), intent(out) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/0,5/)) ) stat= -200

       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(value,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackReal
     
  subroutine PVMIDLunpackSngl(value,info,msg)
    real (r4), intent(out) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    real(r8) :: dble_value

    integer :: details(2), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/0,5/)) ) stat= -200

       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(dble_value,stat)
       value = dble_value
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackSngl
     
  subroutine PVMIDLunpackLogical(value,info,msg)
    logical, intent(out) :: value
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), intValue, stat
    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/=(/0,3/))) stat= -200

       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(intValue,stat)
    end if
    value = intValue /= 0
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackLogical

  subroutine PVMIDLunpackChararr1(values,info,msg)
    character(len=1), intent(out), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(2), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/1,1/)) ) stat= -200 ! rank one byte array

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:1)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackChararr1

  subroutine PVMIDLunpackChararr2(values,info,msg)
    character(len=1), intent(out), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/2,1/)) ) stat= -200 ! rank two byte array

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:2)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackChararr2

  subroutine PVMIDLunpackIntarr1(values,info,msg)
    integer, intent(out), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(2), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/1,3/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:1)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackIntarr1

  subroutine PVMIDLunpackIntarr2(values,info,msg)
    integer, intent(out), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/2,3/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:2)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackIntarr2

  subroutine PVMIDLunpackIntarr3(values,info,msg)
    integer, intent(out), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/3,3/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:3)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackIntarr3
     
  subroutine PVMIDLunpackRealarr1(values,info,msg)
    real (r8), intent(out), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(2), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/1,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:1)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackRealarr1

  subroutine PVMIDLunpackRealarr2(values,info,msg)
    real (r8), intent(out), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/2,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:2)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackRealarr2

  subroutine PVMIDLunpackRealarr3(values,info,msg)
    real (r8), intent(out), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/3,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:3)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(values,stat)
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackRealarr3

  subroutine PVMIDLunpackSnglarr1(values,info,msg)
    real (r4), intent(out), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    real(r8), dimension(size(values)) :: dble_values

    integer :: details(2), sentShape(2), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/1,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:1)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(dble_values,stat)
       values = dble_values
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackSnglarr1

  subroutine PVMIDLunpackSnglarr2(values,info,msg)
    real (r4), intent(out), dimension(:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    real(r8), dimension(size(values,1),size(values,2)) :: dble_values

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/2,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:2)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(dble_values,stat)
       values = dble_values
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackSnglarr2

  subroutine PVMIDLunpackSnglarr3(values,info,msg)
    real (r4), intent(out), dimension(:,:,:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg
    real(r8), dimension(size(values,1),size(values,2),size(values,3)) :: &
      & dble_values

    integer :: details(2), sentShape(3), stat

    ! First unpack noDims and a 5 to indicate double
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/3,5/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:3)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(dble_values,stat)
       values = dble_values
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackSnglarr3

  subroutine PVMIDLunpackLogarr1(values,info,msg)
    logical, intent(out), dimension(:) :: values
    integer, intent(out), optional :: info
    character (LEN=*), intent(in), optional :: msg

    integer :: details(2), sentShape(2), stat
    integer, dimension(size(values)) :: valAsInt

    ! First unpack noDims and a 3 to indicate integer (LONG in IDL of course)
    call pvmf90unpack( details, stat)

    if (stat==0) then 
       if (any(details/= (/1,3/)) ) stat= -200

       ! Now output the dimensions themselves
       if (stat==0) call pvmf90unpack( sentShape,stat)
       if (any(sentShape(1:1)/=shape(values))) stat= -201
       
       ! Now unpack the data itself
       if (stat==0) call pvmf90unpack(valAsInt,stat)
       if (stat==0) values = ( valAsInt == 1 )
    end if
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLunpackLogarr1

  ! ----------------------------------------------------------------------

  ! Now some simple all in one routines to do the initsend, the pack and the
  ! send all at once.  This should make life easier

  subroutine PVMIDLSendString(value,tid,info,msgTag,msg)
    character (LEN=*), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, MYMSGTAG, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid,myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendString

  subroutine PVMIDLSendInteger(value,tid,info,msgTag,msg)
    integer, intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, MYMSGTAG, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value, stat, msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendInteger

  subroutine PVMIDLSendReal(value,tid, info, msgTag, msg)
    real(r8), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendReal

  subroutine PVMIDLSendLogical(value,tid,info,msgTag,msg)
    logical, intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, MYMSGTAG, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value, stat, msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendLogical

  subroutine PVMIDLSendIntArr1(value,tid, info, msgTag, msg)
    integer, dimension(:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendIntArr1

  subroutine PVMIDLSendIntArr2(value,tid, info, msgTag, msg)
    integer, dimension(:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendIntArr2

  subroutine PVMIDLSendIntArr3(value,tid, info, msgTag, msg)
    integer, dimension(:,:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendIntArr3

  subroutine PVMIDLSendRealArr1(value,tid, info, msgTag, msg)
    real(r8), dimension(:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag, stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendRealArr1

  subroutine PVMIDLSendRealArr2(value,tid, info, msgTag, msg)
    real(r8), dimension(:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendRealArr2

  subroutine PVMIDLSendRealArr3(value,tid, info, msgTag, msg)
    real(r8), dimension(:,:,:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag, stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendRealArr3
 

  subroutine PVMIDLSendLogArr1(value,tid, info, msgTag, msg)
    logical, dimension(:), intent(in) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    integer :: bufferID, myMsgTag, stat

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    call PVMFInitSend(PvmDataDefault,bufferID)
    call PVMIDLPack(value,stat,msg)
    if (stat==0) call PVMFSend(tid, myMsgTag,stat)
    call PVMIDLStat ( stat, info, msg )
  end subroutine PVMIDLSendLogArr1

  ! ----------------------------------------------------------------------

  ! Now the same for receive

  subroutine PVMIDLReceiveString(value,tid,info, noBlock, msgTag, msg)
    character (LEN=*), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveString

  subroutine PVMIDLReceiveInteger(value,tid,info, noBlock, msgTag, msg)
    integer, intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveInteger

  subroutine PVMIDLReceiveReal(value,tid,info, noBlock, msgTag, msg)
    real(r8), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

    logical :: useNoBlock = .false.
    integer :: bufferID, myMsgTag

    myMsgTag = IDLMsgTag
    if (present(msgTag)) myMsgTag = msgTag

    if (useNoBlock) then
       call PVMFNrecv(tid, myMsgTag,bufferID)
    else
       call PVMFrecv(tid, myMsgTag,bufferID)
    endif       
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveReal

  subroutine PVMIDLReceiveLogical(value,tid,info, noBlock, msgTag, msg)
    logical, intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveLogical

  subroutine PVMIDLReceiveIntArr1(value,tid,info, noBlock, msgTag, msg)
    integer, dimension(:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveIntArr1

  subroutine PVMIDLReceiveIntArr2(value,tid,info, noBlock, msgTag, msg)
    integer, dimension(:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveIntArr2

  subroutine PVMIDLReceiveIntArr3(value,tid,info, noBlock, msgTag, msg)
    integer, dimension(:,:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveIntArr3

  subroutine PVMIDLReceiveRealArr1(value,tid,info, noBlock, msgTag, msg)
    real(r8), dimension(:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveRealArr1

  subroutine PVMIDLReceiveRealArr2(value,tid,info, noBlock, msgTag, msg)
    real(r8), dimension(:,:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveRealArr2

  subroutine PVMIDLReceiveRealArr3(value,tid,info, noBlock, msgTag, msg)
    real(r8), dimension(:,:,:), intent(out) :: value
    integer, intent(in):: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveRealArr3

  subroutine PVMIDLReceiveLogArr1(value,tid,info, noBlock, msgTag, msg)
    logical, dimension(:), intent(out) :: value
    integer, intent(in) :: tid
    integer, intent(out), optional :: info
    logical, intent(in), optional :: noBlock
    integer, intent(in), optional :: msgTag
    character (LEN=*), intent(in), optional :: msg

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
    call PVMIDLUnpack(value,info,msg)
  end subroutine PVMIDLReceiveLogArr1

! =====     Private procedures     =====================================
  subroutine PVMIDLStat ( stat, info, message )
    integer, intent(in) :: Stat ! Status from PVM operation
    integer, intent(out), optional :: Info ! User's status argument
    character(len=*), intent(in), optional :: Message ! User's error message

    if ( present(info) ) info = stat
    if ( stat /= 0 ) then
      if ( present(message) ) then
        call PVMErrorMessage ( stat, message )
      else if ( .not. present(info) ) then
        call PVMErrorMessage ( stat, &
          & "PVMIDLPack Error occurred and neither INFO nor MESSAGE present" )
      end if
    end if
  end subroutine PVMIDLStat

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module PVMIDL

! $Log$
! Revision 2.12  2005/03/15 23:48:55  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 2.11  2004/10/19 22:59:33  vsnyder
! Add optional 'msg' argument and internal error processing
!
! Revision 2.10  2002/12/04 21:54:58  livesey
! Bug fix in string unpacking
!
! Revision 2.9  2002/10/08 17:43:19  livesey
! Bug fixes
!
! Revision 2.8  2002/10/07 23:22:09  pwagner
! Added (un)packSngl routines
!
