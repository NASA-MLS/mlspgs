! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module QuantityPVM                      ! Send and receive vector quantities using pvm

  ! This module provides functionality for sending and receiving vectors
  ! through a pvm connection.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMFSEND, PVMERRORMESSAGE, &
    & PVMFRECV
  use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
  use String_Table, only: GET_STRING, DISPLAY_STRING
  use VectorsModule, only: VECTORVALUE_T
  use MLSCommon, only: R8
  use Intrinsic, only: LIT_INDICES
  use MLSSignals_m, only: GETSIGNALNAME
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
!---------------------------------------------------------------------------

  public :: PVMSendQuantity, PVMReceiveQuantity

  ! Local parameters
  integer, parameter :: QTYMSGTAG = 200

contains ! ================================== Module procedures ============

  ! ---------------------------------- PVMSendQuantity ---------------------
  subroutine PVMSendQuantity ( Q, tid, justPack, noValues, noMask )
    type (VectorValue_T), intent(in) :: Q ! Quantity to send
    integer, intent(in), optional :: TID ! Task to send it to
    logical, intent(in), optional :: JUSTPACK ! Just pack it into an existing buffer
    logical, intent(in), optional :: NOVALUES ! Don't send the values
    logical, intent(in), optional :: NOMASK ! Don't send the mask

    ! Local variables
    integer :: BUFFERID                 ! From pvm
    integer :: INFO                     ! Flag
    logical :: MYJUSTPACK               ! Copy of justPack
    logical :: MYNOVALUES               ! Copy of noValues
    logical :: MYNOMASK                 ! Copy of noMask

    character(len=132) :: WORD          ! Result of get_string etc.

    ! Executable code
    myJustPack = .false.
    myNoValues = .false.
    myNoMask = .false.
    if ( present(justPack) ) myJustPack = justPack
    if ( present(noValues) ) myNoValues = noValues
    if ( present(noMask) ) myNoMask = noMask

    ! Now we simply pack the quantity up and send it down the pvm spigot
    if (.not. myJustPack) call PVMFInitSend ( PvmDataDefault, bufferID )
    
    call PVMIDLPack ( (/ q%template%noInstances, &
      & q%template%noSurfs, q%template%noChans, q%template%instanceLen /), &
      & info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity dimensions." )

    call PVMIDLPack ( (/ q%template%coherent, &
      & q%template%stacked, q%template%regular, q%template%minorFrame, &
      & q%template%logBasis /), info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity flags" )

    call PVMIDLPack ( (/ q%template%noInstancesLowerOverlap, &
      & q%template%noInstancesUpperOverlap, q%template%sideband, &
      & q%template%instrumentModule, q%template%radiometer, &
      & q%template%quantityType, q%template%unit, q%template%frequencyCoordinate, &
      & q%template%molecule, q%template%verticalCoordinate,&
      & q%template%signal, q%template%name /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing misc quantity stuff" )

    ! Now pack some strings

    call Get_String( lit_indices(q%template%quantityType), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantityType" )

    call Get_String( lit_indices(q%template%verticalCoordinate), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing verticalCoordinate" )

    call Get_String( lit_indices(q%template%unit), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing unit" )

    call Get_String( lit_indices(q%template%frequencyCoordinate), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing frequencyCoordinate" )

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

    if ( associated ( q%template%mafIndex ) ) then
      call PVMIDLPack ( (/.true./), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing flag" )
      call PVMIDLPack ( q%template%mafIndex, info )
    else
      call PVMIDLPack ( (/ .false. /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing mafIndex/flag" )

    if ( associated ( q%template%mafCounter ) ) then
      call PVMIDLPack ( (/.true./), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing flag" )
      call PVMIDLPack ( q%template%mafCounter, info )
    else
      call PVMIDLPack ( (/ .false. /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing mafCounter/flag" )

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
  subroutine PVMReceiveQuantity ( QT, values, tid, mask, justUnpack )
    type (QuantityTemplate_T), intent(out) :: QT ! Template for quantity
    real (r8), dimension(:,:), pointer :: VALUES ! Values for quantity
    integer, intent(in), optional :: TID ! Task to get it from
    integer, dimension(:,:), optional, pointer :: MASK ! Mask
    logical, intent(in), optional :: JUSTUNPACK ! Just unpack from existing buffer

    ! Local variables
    integer :: BUFFERID                 ! From pvm
    integer :: INFO                     ! Flag
    integer :: I4(4)                    ! Unpacked stuff
    integer :: I12(12)                  ! Unpacked stuff
    logical :: L2(2)                    ! Unpacked stuff
    logical :: L5(5)                    ! Unpacked stuff
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

    call PVMIDLUnPack ( l5, info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking quantity flags" )

    call SetupNewQuantityTemplate ( qt, &
      & noInstances  = i4(1), &
      & noSurfs      = i4(2), &
      & noChans      = i4(3), &
      & coherent     = l5(1), &
      & stacked      = l5(2), &
      & regular      = l5(3), &
      & instanceLen  = i4(4), &
      & minorFrame   = l5(4) )
    qt%logBasis = l5(5)

    call PVMIDLUnPack ( i12, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking misc quantity stuff" )
    qt%noInstancesLowerOverlap = i12(1)
    qt%noInstancesUpperOverlap = i12(2)
    qt%sideband                = i12(3)
    qt%instrumentModule        = i12(4)
    qt%radiometer              = i12(5)
    qt%quantityType            = i12(6)
    qt%unit                    = i12(7)
    qt%frequencyCoordinate     = i12(8)
    qt%molecule                = i12(9)
    qt%verticalCoordinate      = i12(10)
    qt%signal                  = i12(11)
    qt%name                    = i12(12)

    ! Now unpack some strings
    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking quantityType" )
    ! Just ignore it, we got it as an integer already
    
    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking verticalCoordinate" )
    ! Just ignore it, we got it as an integer already

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking unit" )
    ! Just ignore it, we got it as an integer already

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking frequencyCoordiante" )
    ! Just ignore it, we got it as an integer already

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking signal" )
    ! Just ignore it, we got it as an integer already

    call PVMIDLUnpack ( word, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking molecule" )
    ! Juest ignore it, we got it as an integer already

    ! Now unpack the arrays

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


    call PVMIDLUnpack ( flag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking flag" )
    if ( flag(1) ) then
      call PVMIDLUnpack ( qt%mafIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking mafIndex" )
    end if

    call PVMIDLUnpack ( flag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking flag" )
    if ( flag(1) ) then
      call PVMIDLUnpack ( qt%mafCounter, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking mafCounter" )
    end if

    call PVMIDLUnpack ( flag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking flag" )
    if ( flag(1) ) then
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
    nullify ( values )
    call Allocate_Test ( values, qt%instanceLen, qt%noInstances, &
      & 'values', ModuleName )

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
    if ( .not. noMask .and. present(mask) ) then
      call CreateMaskArray ( mask, values )
      call PVMIDLUnpack ( mask, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "unpacking mask" )
    end if
      
  end subroutine PVMReceiveQuantity

end module QuantityPVM


