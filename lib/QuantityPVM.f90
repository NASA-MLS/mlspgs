! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

module QuantityPVM                      ! Send and receive vector quantities using pvm
  ! This module provides functionality for sending and receiving vectors
  ! through a pvm connection.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFSEND, &
    & PVMFRECV
  use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
  use MorePVM, only: PVMPackLitIndex, PVMPackStringIndex, &
    & PVMUnpackLitIndex, PVMUnpackStringIndex
  use String_Table, only: GET_STRING, DISPLAY_STRING
  use VectorsModule, only: VECTORVALUE_T
  use MLSCommon, only: R8
  use Intrinsic, only: LIT_INDICES
  use MLSSignals_m, only: GETSIGNALNAME
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, PVMERRORMESSAGE
  use QuantityTemplates, only: QUANTITYTEMPLATE_T, SETUPNEWQUANTITYTEMPLATE
  use Dump_0, only: DUMP
  use VectorsModule, only: CREATEMASKARRAY

  implicit none
  private

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  public :: PVMSendQuantity, PVMReceiveQuantity

  ! Local parameters
  integer, parameter :: QTYMSGTAG = 200

contains ! ================================== Module procedures ============

  ! ---------------------------------- PVMSendQuantity ---------------------
  subroutine PVMSendQuantity ( Q, tid, justPack, noValues, noMask, skipMIFGeolocation )
    type (VectorValue_T), intent(in) :: Q ! Quantity to send
    integer, intent(in), optional :: TID ! Task to send it to
    logical, intent(in), optional :: JUSTPACK ! Just pack it into an existing buffer
    logical, intent(in), optional :: NOVALUES ! Don't send the values
    logical, intent(in), optional :: NOMASK ! Don't send the mask
    logical, intent(in), optional :: SKIPMIFGEOLOCATION ! Skip geolocation if minor frame

    ! Local variables
    integer :: BUFFERID                 ! From pvm
    integer :: INFO                     ! Flag
    logical :: MYJUSTPACK               ! Copy of justPack
    logical :: MYNOVALUES               ! Copy of noValues
    logical :: MYNOMASK                 ! Copy of noMask
    logical :: MYSKIPMIFGEOLOCATION     ! If set, skip geolocation

    character(len=132) :: WORD          ! Result of get_string etc.

    ! Executable code
    myJustPack = .false.
    myNoValues = .false.
    myNoMask = .false.
    mySkipMIFGeolocation = .false.
    if ( present(justPack) ) myJustPack = justPack
    if ( present(noValues) ) myNoValues = noValues
    if ( present(noMask) ) myNoMask = noMask
    if ( present(skipMIFGeolocation) ) mySkipMIFGeolocation = skipMIFGeolocation

    ! Now we simply pack the quantity up and send it down the pvm spigot
    if (.not. myJustPack) call PVMFInitSend ( PvmDataDefault, bufferID )
    
    call PVMIDLPack ( (/ q%template%noInstances, &
      & q%template%noSurfs, q%template%noChans, q%template%instanceLen /), &
      & info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity dimensions." )

    call PVMIDLPack ( (/ q%template%coherent, &
      & q%template%stacked, q%template%regular, q%template%minorFrame, &
      & q%template%majorFrame, q%template%logBasis /), info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity flags" )

    call PVMIDLPack ( (/ q%template%noInstancesLowerOverlap, &
      & q%template%noInstancesUpperOverlap, q%template%sideband, &
      & q%template%instrumentModule, q%template%radiometer, &
      & q%template%instanceOffset, q%template%grandTotalInstances /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing misc quantity stuff" )

    ! Now pack some strings

    call PVMPackLitIndex ( q%template%quantityType, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantityType" )
    call PVMPackLitIndex ( q%template%verticalCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing verticalCoordinate" )
    call PVMPackLitIndex ( q%template%frequencyCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing frequencyCoordinate" )
    call PVMPackLitIndex ( q%template%reflector, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing reflector" )
    call PVMPackLitIndex ( q%template%unit, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing unit" )
    call PVMPackLitIndex ( q%template%molecule, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing molecule" )
    call PVMPackStringIndex ( q%template%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing name" )

    ! Pack signal as a string
    if ( q%template%signal /= 0 ) then
      call GetSignalName( q%template%signal, &
        & word, sideband=q%template%sideband )
    else
      word = ''
    endif
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing signal" )

    if ( q%template%molecule >= lbound(lit_indices,1) .and. &
      &  q%template%molecule <= ubound(lit_indices,1) ) then
      call Get_String( lit_indices(q%template%molecule), word, noError=.true. )
    else
      word = ''
    endif
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing molecule" )

    ! Now pack the arrays
    if ( .not. (q%template%minorFrame .or. q%template%majorFrame) &
      & .or. .not. mySkipMIFGeolocation ) then
      call PVMIDLPack ( q%template%surfs, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing surfs" )
      
      call PVMIDLPack ( q%template%phi, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing phi" )
      
      call PVMIDLPack ( q%template%geodLat, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing geodLat" )
      
      call PVMIDLPack ( q%template%lon, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing lon" )
      
      call PVMIDLPack ( q%template%time, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing time" )
      
      call PVMIDLPack ( q%template%solarTime, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing solarTime" )
      
      call PVMIDLPack ( q%template%solarZenith, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing solarZenith" )
      
      call PVMIDLPack ( q%template%losAngle, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing losAngle" )
    end if

    if ( associated ( q%template%frequencies ) ) then
      call PVMIDLPack ( (/.true./), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing flag" )
      call PVMIDLPack ( q%template%frequencies, info )
    else
      call PVMIDLPack ( (/ .false. /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing frequencies/flag" )

    ! Finally the two arrays for irregular quantities

    if ( .not. q%template%regular ) then
      call PVMIDLPack ( q%template%surfIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing surfIndex" )

      call PVMIDLPack ( q%template%chanIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing chanIndex" )
    end if

    ! Now pack the noValues, noMask flag
    myNoMask = myNoMask .or. ( .not. associated(q%mask) )
    call PVMIDLPAck ( (/ myNoValues, myNoMask /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing values/mask flags" )    

    ! Now we're going to send this to the snooper.

    if (.not. myJustPack) then
      call PVMFSend ( tid, QtyMsgTag, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector template" )
      
      ! Now we're going to send the values in a separate message
      if ( .not. all ( (/myNoValues, myNoMask /) ) ) &
        & call PVMFInitSend ( PVMDataDefault, bufferID)
    end if

    ! Pack the values
    if ( .not. myNoValues ) then 
      call PVMIDLPack ( q%values, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending values" )
    endif

    ! Skip the mask for the moment.
    if ( .not. myNoMask ) then
      call PVMIDLPack ( q%mask, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending mask" )
    endif

    ! Send this buffer
    if (.not. myJustPack .and. .not. all ((/myNoValues, myNoMask/)) ) then
      call PVMFSend ( tid, QtyMsgTag, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector values" )
    end if
      
  end subroutine PVMSendQuantity

  ! ---------------------------------- PVMReceiveQuantity ---------------------
  subroutine PVMReceiveQuantity ( QT, values, tid, mask, justUnpack, mifGeolocation )
    type (QuantityTemplate_T), intent(out) :: QT ! Template for quantity
    ! It's not inout, because then setupNewQuantityTemplate would deallocate
    ! the pointer components.  But the actual argument is put into a database
    ! using a shallow copy, so cleaning it up would clobber a database item.
    real (r8), dimension(:,:), optional, pointer :: VALUES ! Values for quantity
    integer, intent(in), optional :: TID ! Task to get it from
    character, dimension(:,:), optional, pointer :: MASK ! Mask
    logical, intent(in), optional :: JUSTUNPACK ! Just unpack from existing buffer
    type (QuantityTemplate_T), dimension(:), intent(in), optional :: MIFGEOLOCATION

    ! Local variables
    integer :: BUFFERID                 ! From pvm
    integer :: INFO                     ! Flag
    integer :: I4(4)                    ! Unpacked stuff
    integer :: I7(7)                    ! Unpacked stuff
    logical :: L2(2)                    ! Unpacked stuff
    logical :: L6(6)                    ! Unpacked stuff
    logical :: FLAG(1)                  ! To unpack
    character(len=132) :: WORD          ! Result of get_string etc.
    logical :: MYJUSTUNPACK             ! Copy of justUnPack
    logical :: NOVALUES                 ! No values sent
    logical :: NOMASK                   ! No mask sent

    ! Executable code

    myJustUnpack = .false.
    if ( present (justUnpack) ) myJustUnPack = justUnpack

    if ( .not. myJustUnpack ) then
      ! Get buffer, we'll wait for it, assume the calling code knows it's coming.
      call PVMFrecv ( tid, QtyMsgTag, bufferID )
    end if

    ! Now we unpack the information

    ! Get dimensions
    call PVMIDLUnpack( i4, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking quantity dimensions." )

    call PVMIDLUnPack ( l6, info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking quantity flags" )

    call SetupNewQuantityTemplate ( qt, &
      & noInstances  = i4(1), &
      & noSurfs      = i4(2), &
      & noChans      = i4(3), &
      & coherent     = l6(1), &
      & stacked      = l6(2), &
      & regular      = l6(3), &
      & instanceLen  = i4(4), &
      & minorFrame   = l6(4), &
      & majorFrame   = l6(5) )
    qt%logBasis = l6(6)

    call PVMIDLUnPack ( i7, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking misc quantity stuff" )
    qt%noInstancesLowerOverlap = i7(1)
    qt%noInstancesUpperOverlap = i7(2)
    qt%sideband                = i7(3)
    qt%instrumentModule        = i7(4)
    qt%radiometer              = i7(5)
    qt%instanceOffset          = i7(6)
    qt%grandTotalInstances     = i7(7)

    ! Now unpack literals etc.
    call PVMUnpackLitIndex ( qt%quantityType, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking quantityType" )
    call PVMUnpackLitIndex ( qt%verticalCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking verticalCoordinate" )
    call PVMUnpackLitIndex ( qt%frequencyCoordinate, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking frequencyCoordiante" )
    call PVMUnpackLitIndex ( qt%unit, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking unit" )
    call PVMUnpackLitIndex ( qt%reflector, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking reflector" )
    call PVMUnpackLitIndex ( qt%molecule, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking molecule" )
    call PVMUnpackStringIndex ( qt%name, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking name" )

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking signal" )
    ! Just ignore it, we got it as an integer already

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking molecule" )
    ! Juest ignore it, we got it as an integer already

    ! Now unpack the arrays
    if ( .not. ( qt%minorFrame .or. qt%majorFrame ) .or. &
      & .not. present ( mifGeolocation ) ) then

      call PVMIDLUnpack ( qt%surfs, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking surfs" )

      call PVMIDLUnpack ( qt%phi, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking phi" )

      call PVMIDLUnpack ( qt%geodLat, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking geodLat" )

      call PVMIDLUnpack ( qt%lon, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking lon" )

      call PVMIDLUnpack ( qt%time, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking time" )

      call PVMIDLUnpack ( qt%solarTime, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking solarTime" )

      call PVMIDLUnpack ( qt%solarZenith, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking solarZenith" )

      call PVMIDLUnpack ( qt%losAngle, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking losAngle" )

    else
      ! If it's minor frame and we've got mif geolocation information just point to that.
      if ( qt%minorFrame ) then
        qt%surfs => mifGeolocation(qt%instrumentModule)%surfs
        qt%phi => mifGeolocation(qt%instrumentModule)%phi
        qt%geodLat => mifGeolocation(qt%instrumentModule)%geodLat
        qt%lon => mifGeolocation(qt%instrumentModule)%lon
        qt%time => mifGeolocation(qt%instrumentModule)%time
        qt%solarTime => mifGeolocation(qt%instrumentModule)%solarTime
        qt%solarZenith => mifGeolocation(qt%instrumentModule)%solarZenith
      end if
      if ( qt%majorFrame ) then
        nullify ( qt%surfs )
        qt%phi => mifGeolocation(qt%instrumentModule)%phi(1:1,:)
        qt%geodLat => mifGeolocation(qt%instrumentModule)%geodLat(1:1,:)
        qt%lon => mifGeolocation(qt%instrumentModule)%lon(1:1,:)
        qt%time => mifGeolocation(qt%instrumentModule)%time(1:1,:)
        qt%solarTime => mifGeolocation(qt%instrumentModule)%solarTime(1:1,:)
        qt%solarZenith => mifGeolocation(qt%instrumentModule)%solarZenith(1:1,:)
      end if
    end if

    call PVMIDLUnpack ( flag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking flag" )
    if ( flag(1) ) then
      nullify ( qt%frequencies )
      call Allocate_test ( qt%frequencies, qt%noChans, &
        & 'qt%frequencies', ModuleName )
      call PVMIDLUnpack ( qt%frequencies, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking frequencies" )
    end if

    ! Finally the two arrays for irregular quantities
    if ( .not. qt%regular ) then
      call PVMIDLUnpack ( qt%surfIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking surfIndex" )

      call PVMIDLUnpack ( qt%chanIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking chanIndex" )
    end if

    ! Now the value/mask flags
    call PVMIDLUnPack ( l2, info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking value/mask flags" )
    noValues = l2(1)
    noMask = l2(2)

    ! Setup the values
    if ( .not. noValues .and. .not. present ( values ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Values sent but no place to put them' )

    if ( present ( values ) ) then
      nullify ( values )
      call Allocate_Test ( values, qt%instanceLen, qt%noInstances, &
        & 'values', ModuleName )
    endif

    if ( .not. myJustUnpack .and. .not. all ((/noValues,noMask/)) ) then
      ! Now we're going to receive the values in a separate message
      call PVMFrecv ( tid, QtyMsgTag, bufferID )
    end if

    ! Unpack the values
    if ( .not. noValues ) then
      call PVMIDLUnpack ( values, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking values" )
    endif

    ! Skip the mask for the moment.
    if ( present(mask) ) then
      nullify ( mask )
      if ( .not. noMask ) then
        call CreateMaskArray ( mask, values )
        call PVMIDLUnpack ( mask, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking mask" )
      endif
    else
      if ( .not. noMask ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Mask sent but no place to put it' )
    end if

  end subroutine PVMReceiveQuantity

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module QuantityPVM

! $Log$
! Revision 2.19  2005/03/15 23:48:55  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule
!
! Revision 2.18  2003/07/07 20:28:53  livesey
! Removed a print statement
!
! Revision 2.17  2003/07/07 20:22:00  livesey
! Transfers additional stuff from quantity template
!
! Revision 2.16  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.15  2003/05/13 04:46:55  livesey
! Changed some integer packing to strings
!
! Revision 2.14  2002/11/27 00:24:01  livesey
! Better handling of major frame quantities
!
! Revision 2.13  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/10/06 01:08:53  livesey
! Put in the skipMIFGeolocation functionality
!
! Revision 2.11  2002/10/05 00:42:52  livesey
! Modified to use stuff from MorePVM
!
! Revision 2.10  2002/08/16 21:41:13  livesey
! Bug fix in frequency transmission
!
! Revision 2.9  2002/07/01 23:50:46  vsnyder
! Add an important comment about memory leakage
!
! Revision 2.8  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.7  2002/02/01 23:50:36  livesey
! Added CVS log information
!
