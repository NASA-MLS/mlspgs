module IDLCFM2_m
    use CFM
    use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
    use PVM, only: PVMFRECV, PVMFINITSEND, PVMDATADEFAULT
    use MorePVM, only: PVMUnpackStringIndex, PVMUnpackLitIndex
    use MLSMessageModule, only: PVMERRORMESSAGE
    use Allocate_Deallocate, only: ALLOCATE_TEST, Deallocate_test

    implicit none
    private
    public :: ICFMReceiveQuantity, ICFMReceiveVector, QTYMSGTAG, ICFMSendVector
!---------------------------- RCS Ident Info -------------------------------
    character(len=*), private, parameter :: ModuleName= &
        "$RCSfile$"
    private :: not_used_here
!---------------------------------------------------------------------------

    ! Local parameters
    integer, parameter :: QTYMSGTAG = 200

    contains
!--------------------------- end bloc --------------------------------------
    logical function not_used_here()
        character (len=*), parameter :: ModuleName= &
            "$RCSfile$"
        character (len=*), parameter :: IdParm = &
            "$Id$"
        character (len=len(idParm)) :: Id = idParm
        not_used_here = (id(1:1) == ModuleName(1:1))
        print *, Id ! .mod files sometimes change if PRINT is added
    end function not_used_here
!---------------------------------------------------------------------------

    subroutine ICFMSendQuantity (qty, info, tid, callsend)
        type(VectorValue_T), intent(in) :: qty
        integer, intent(out) :: info
        integer, intent(in), optional :: tid
        logical, intent(in), optional :: callsend

        logical :: mycallsend = .true.
        integer :: bufid
        logical :: l28(28)
        type(QuantityTemplate_T) :: template

        integer, parameter :: P_NAME = 1
        integer, parameter :: P_TYPE = 2
        integer, parameter :: P_OFFSET = 3
        integer, parameter :: P_COHERENCE = 4
        integer, parameter :: P_LOGBASIS = 5
        integer, parameter :: P_MINVALUE = 6
        integer, parameter :: P_BADVALUE = 7
        integer, parameter :: P_MOLECULE = 8
        integer, parameter :: P_MODULE = 9
        integer, parameter :: P_SIGNAL = 10
        integer, parameter :: P_RADIOMETER = 11
        integer, parameter :: P_FREQUENCIES = 12
        integer, parameter :: P_FCOORD = 13
        integer, parameter :: P_NOCHANS = 14
        integer, parameter :: P_SURFS = 15
        integer, parameter :: P_VCOORD = 16
        integer, parameter :: P_NOSURFS = 17
        integer, parameter :: P_PHI = 18
        integer, parameter :: P_GEODLAT = 19
        integer, parameter :: P_LONGITUDE = 20
        integer, parameter :: P_LOSANGLE = 21
        integer, parameter :: P_SOLARZENITH = 22
        integer, parameter :: P_SOLARTIME = 23
        integer, parameter :: P_TIME = 24
        integer, parameter :: P_NOPROFS = 25
        integer, parameter :: P_INSTANCELEN = 26
        integer, parameter :: P_VALUE = 27
        integer, parameter :: P_MASK = 28

        if (present(callsend) .and. .not. callsend) mycallsend = .false.

        if (mycallsend .and. .not. present(tid)) call MLSMessage (MLSMSG_Error, moduleName, "Missing 'tid'")

        if (mycallsend) call PVMFInitSend(PVMDataDefault, bufid)

        template = qty%template

        l28(P_NAME) = template%name == 0
        l28(P_TYPE) = template%quantityType == 0
        l28(P_OFFSET) = .true.
        l28(P_COHERENCE) = .true.
        l28(P_LOGBASIS) = .true.
        l28(P_MINVALUE) = template%logBasis !because if logBasis is false, no use for minValue
        l28(P_BADVALUE) = .true.
        l28(P_MOLECULE) = template%molecule == 0
        l28(P_MODULE) = template%instrumentModule == 0
        l28(P_SIGNAL) = template%signal == 0
        l28(P_RADIOMETER) = template%radiometer == 0
        l28(P_FREQUENCIES) = associated(template%frequencies) .and. size(template%frequencies) .gt. 0
        l28(P_FCOORD) = l28(P_FREQUENCIES)
        l28(P_NOCHANS) = template%noChans == 0
        l28(P_NOSURFS) = template%noSurfs == 0
        l28(P_SURFS) = associated(template%surfs) .and. size(template%surfs) .gt. 0
        l28(P_VCOORD) = l28(P_SURFS)
        l28(P_PHI) = associated(template%phi) .and. size(template%phi) .gt. 0
        l28(P_GEODLAT) = associated(template%geodlat) .and. size(template%geodlat) .gt. 0
        l28(P_LONGITUDE) = associated(template%lon) .and. size(template%lon) .gt. 0
        l28(P_LOSANGLE) = associated(template%losAngle) .and. size(template%losangle) .gt. 0
        l28(P_SOLARZENITH) = associated(template%solarZenith) .and. size(template%solarzenith) .gt. 0
        l28(P_SOLARTIME) = associated(template%solarTime) .and. size(template%solartime) .gt. 0
        l28(P_TIME) = associated(template%time) .and. size(template%time) .gt. 0
        ! this variable reflect the size of multiple arrays, so just send it for simplicity
        l28(P_NOPROFS) = .true.
        l28(P_INSTANCELEN) = .true.
        l28(P_VALUE) = associated(qty%values) .and. size(qty%values) .gt. 0
        l28(P_MASK) = associated(qty%mask) .and. size(qty%mask) .gt. 0

        call PVMIDLPack(l28, info)
        print *, "packing l28 ", l28

        if (mycallsend) then
            call PVMFSend ( tid, QtyMsgTag, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector value" )
        endif
    end subroutine

    subroutine ICFMReceiveQuantity ( QT, values, mask, tid, callrecv)
        use Allocate_Deallocate, only: ALLOCATE_TEST
        use MLSCommon, only: R8
        use QuantityTemplates, only: QUANTITYTEMPLATE_T, SETUPNEWQUANTITYTEMPLATE
        use Intrinsic, only: l_none
        use VectorsModule, only: CREATEMASKARRAY
        use ConstructQuantityTemplates, only: noProperties, propertyTable, p_majorFrame, &
                                              p_minorFrame, unitstable
        use parse_signal_m, only: parse_signal
        use MLSSignals_m, only: GetRadiometerIndex, GetModuleFromRadiometer, &
                                GetModuleIndex, GetRadiometerFromSignal, &
                                GetModuleFromSignal

        type (QuantityTemplate_T), intent(out) :: QT ! Template for quantity
        ! It's not inout, because then setupNewQuantityTemplate would deallocate
        ! the pointer components.  But the actual argument is put into a database
        ! using a shallow copy, so cleaning it up would clobber a database item.
        real (r8), dimension(:,:), pointer :: VALUES ! Values for quantity
        character, dimension(:,:), pointer :: MASK ! Mask
        integer, intent(in), optional :: TID ! Task to get it from
        logical, optional, intent(in) :: callrecv !true if this subroutine should call PVMFRecv, default is false

        integer, parameter :: P_NAME = 1
        integer, parameter :: P_TYPE = 2
        integer, parameter :: P_OFFSET = 3
        integer, parameter :: P_COHERENCE = 4
        integer, parameter :: P_LOGBASIS = 5
        integer, parameter :: P_MINVALUE = 6
        integer, parameter :: P_BADVALUE = 7
        integer, parameter :: P_MOLECULE = 8
        integer, parameter :: P_MODULE = 9
        integer, parameter :: P_SIGNAL = 10
        integer, parameter :: P_RADIOMETER = 11
        integer, parameter :: P_FREQUENCIES = 12
        integer, parameter :: P_FCOORD = 13
        integer, parameter :: P_NOCHANS = 14
        integer, parameter :: P_SURFS = 15
        integer, parameter :: P_VCOORD = 16
        integer, parameter :: P_NOSURFS = 17
        integer, parameter :: P_PHI = 18
        integer, parameter :: P_GEODLAT = 19
        integer, parameter :: P_LONGITUDE = 20
        integer, parameter :: P_LOSANGLE = 21
        integer, parameter :: P_SOLARZENITH = 22
        integer, parameter :: P_SOLARTIME = 23
        integer, parameter :: P_TIME = 24
        integer, parameter :: P_NOPROFS = 25
        integer, parameter :: P_INSTANCELEN = 26
        integer, parameter :: P_VALUE = 27
        integer, parameter :: P_MASK = 28

        ! Local variables
        integer :: BUFFERID                 ! From pvm
        integer :: INFO                     ! Flag
        logical :: l28(28)
        logical, dimension(noProperties) :: PROPERTIES ! Properties for this quantity type
        character(len=32) :: signalString
        integer, dimension(:), pointer :: SignalInds ! From parse signal
        integer :: sideband
        character(len=16) :: radiometerString
        integer, dimension(2) :: hshape ! shape of hgrid-related fields

        ! Executable code
        hshape = 0

        ! Get buffer, we'll wait for it, assume the calling code knows it's coming.
        if (present(callrecv) .and. callrecv) call PVMFrecv ( tid, QtyMsgTag, bufferID )

        call PVMIDLUnpack(l28, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking l28." )
            return
        endif
!        print *, "l28 ", l28

        if (l28(P_NAME)) then
            call PVMUnpackStringIndex ( qt%name, info)
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking name" )
                return
            endif
!            print *, "name " , qt%name
        endif

        if (l28(P_TYPE)) then
            call PVMUnpackLitIndex ( qt%quantityType, info )
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking quantityType" )
                return
            endif
!            print *, "type " , qt%quantityType
            qt%unit = unitsTable(qt%quantityType)
!            print *, "unit ", qt%unit
            properties = propertyTable(:, qt%quantityType)
            qt%majorFrame = properties(p_majorFrame)
!            print *, "majorFrame ", qt%majorFrame
            qt%minorFrame = properties(p_minorFrame)
!            print *, "minorFrame ", qt%minorFrame
        endif

        if (l28(P_OFFSET)) then
            call PVMIDLUnpack(qt%instanceOffset, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instanceOffset" )
                return
            endif
!            print *, "instanceOffset ", qt%instanceOffset
        endif

        if (l28(P_COHERENCE)) then
            call PVMIDLUnpack(qt%coherent, info)
            if (info /= 0) then
                call PVMErrorMessage(info, "unpacking coherence")
                return
            endif
!            print *, "coherence ", qt%coherent
        endif

        if (l28(P_LOGBASIS)) then
            call PVMIDLUnpack(qt%logbasis, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking logbasis" )
                return
            endif
!            print *, "logbasis ", qt%logbasis
        endif

        if (l28(P_MINVALUE)) then
            call PVMIDLUnpack(qt%minvalue, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking minvalue" )
                return
            endif
!            print *, "minvalue ", qt%minvalue
        endif

        if (l28(P_BADVALUE)) then
            call PVMIDLUnpack(qt%badvalue, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking badvalue" )
                return
            endif
!            print *, "badvalue ", qt%badvalue
        endif

        if (l28(P_MOLECULE)) then
            call PVMUnpackLitIndex(qt%molecule, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking molecule" )
                return
            endif
!            print *, "molecule ", qt%molecule
        endif

        if (l28(P_MODULE)) then
            call PVMIDLUnpack(signalString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instrumentModule" )
                return
            endif
            call GetModuleIndex(signalString, qt%instrumentModule)
!            print *, "instrumentModule ", qt%instrumentModule
        endif

        if (l28(P_SIGNAL)) then
            call PVMIDLUnpack(signalString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking signalString" )
                return
            endif
!            print *, "signalString ", signalString
            nullify(signalInds)
            call parse_Signal ( signalString, signalInds, sideband=sideband)
            qt%signal = signalInds(1)
            qt%sideband = sideband
            call deallocate_test ( signalInds, 'signalInds', ModuleName )
!            print *, "signal ", qt%signal
            qt%radiometer = GetRadiometerFromSignal(qt%signal)
            qt%instrumentModule = GetModuleFromSignal(qt%signal)
        endif

        if (l28(P_RADIOMETER)) then
            call PVMIDLUnpack(radiometerString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking radiometerString" )
                return
            endif
!            print *, "radiometerString ", radiometerString
            call GetRadiometerIndex(radiometerString, qt%radiometer)
            ! Every radiometer pairs with only one instrument module
            qt%instrumentModule = GetModuleFromRadiometer(qt%radiometer)
            ! If the user happen to pass in the wrong instrumentModule,
            ! silently correct it, for simplicity
!            print *, "radiometer ", qt%radiometer
!            print *, "module ", qt%instrumentModule
        endif

        if (l28(P_NOCHANS)) then
            call PVMIDLUnpack(qt%noChans, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noChans" )
                return
            endif
!            print *, "noChans ", qt%noChans
        endif

        if (l28(P_FCOORD)) then
            call PVMUnpackLitIndex(qt%frequencyCoordinate, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking frequencyCoordinate" )
                return
            endif
!            print *, "frequencyCoordinate ", qt%frequencyCoordinate
        endif

        if (l28(P_FREQUENCIES)) then
            call allocate_test(qt%frequencies, qt%noChans, "qt%frequencies", moduleName)
            call PVMIDLUnpack(qt%frequencies, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking frequencies" )
                return
            endif
!            print *, "frequencies ", qt%frequencies
        endif

        if (l28(P_NOSURFS)) then
            call PVMIDLUnpack(qt%noSurfs, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noSurfs" )
                return
            endif
!            print *, "noSurfs ", qt%noSurfs
        endif

        if (l28(P_NOPROFS)) then
            call PVMIDLUnpack(qt%noInstances, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noInstances" )
                return
            endif
!            print *, "noProfs ", qt%noInstances
        endif

        if (l28(P_PHI)) then
            if (qt%coherent) then
                call allocate_test(qt%phi, 1, qt%noInstances, "qt%phi", moduleName)
            else
                call allocate_test(qt%phi, qt%nosurfs, qt%noInstances, "qt%phi", moduleName)
            endif
            call PVMIDLUnpack(qt%phi, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking phi" )
                return
            endif
!            print *, "phi ", qt%phi
            hshape = shape(qt%phi)
        endif

        if (l28(P_GEODLAT)) then
            if (qt%coherent) then
                call allocate_test(qt%geodlat, 1, qt%noInstances, "qt%geodlat", moduleName)
            else
                call allocate_test(qt%geodlat, qt%nosurfs, qt%noInstances, "qt%geodlat", moduleName)
            endif
            call PVMIDLUnpack(qt%geodlat, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking geodlat" )
                return
            endif
!            print *, "geodlat ", qt%geodlat
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%geodlat))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "geodlat's shape is not the same as others' shape")
            endif
            hshape = shape(qt%geodlat)
        endif

        if (l28(P_LONGITUDE)) then
            if (qt%coherent) then
                call allocate_test(qt%lon, 1, qt%noInstances, "qt%lon", moduleName)
            else
                call allocate_test(qt%lon, qt%nosurfs, qt%noInstances, "qt%lon", moduleName)
            endif
            call PVMIDLUnpack(qt%lon, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking lon" )
                return
            endif
!            print *, "lon ", qt%lon
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%lon))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "lon's shape is not the same as others' shape")
            endif
            hshape = shape(qt%lon)
        endif

        if (l28(P_LOSANGLE)) then
            if (qt%coherent) then
                call allocate_test(qt%losAngle, 1, qt%noInstances, "qt%losAngle", moduleName)
            else
                call allocate_test(qt%losAngle, qt%nosurfs, qt%noInstances, "qt%losAngle", moduleName)
            endif
            call PVMIDLUnpack(qt%losAngle, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking losAngle" )
                return
            endif
!            print *, "losAngle ", qt%losAngle
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%losAngle))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "losAngle's shape is not the same as others' shape")
            endif
            hshape = shape(qt%losAngle)
        endif

        if (l28(P_SOLARZENITH)) then
            if (qt%coherent) then
                call allocate_test(qt%solarZenith, 1, qt%noInstances, "qt%solarZenith", moduleName)
            else
                call allocate_test(qt%solarZenith, qt%nosurfs, qt%noInstances, "qt%solarZenith", moduleName)
            endif
            call PVMIDLUnpack(qt%solarZenith, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking solarZenith" )
                return
            endif
!            print *, "solarZenith ", qt%solarZenith
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%solarZenith))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "solarZenith's shape is not the same as others' shape")
            endif
            hshape = shape(qt%solarZenith)
        endif

        if (l28(P_SOLARTIME)) then
            if (qt%coherent) then
                call allocate_test(qt%solartime, 1, qt%noInstances, "qt%solartime", moduleName)
            else
                call allocate_test(qt%solartime, qt%nosurfs, qt%noInstances, "qt%solartime", moduleName)
            endif
            call PVMIDLUnpack(qt%solartime, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking solartime" )
                return
            endif
!            print *, "solartime ", qt%solartime
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%solartime))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "solarTime's shape is not the same as others' shape")
            endif
            hshape = shape(qt%solartime)
        endif

        if (l28(P_TIME)) then
            if (qt%coherent) then
                call allocate_test(qt%time, 1, qt%noInstances, "qt%time", moduleName)
            else
                call allocate_test(qt%time, qt%nosurfs, qt%noInstances, "qt%time", moduleName)
            endif
            call PVMIDLUnpack(qt%time, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking time" )
                return
            endif
!            print *, "time ", qt%time
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%time))) then
                call MLSMessage(MLSMSG_Warning, moduleName, "time's shape is not the same as others' shape")
            endif
            hshape = shape(qt%time)
        endif

        if (l28(P_VCOORD)) then
            call PVMUnpackLitIndex(qt%verticalCoordinate, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking verticalCoordinate" )
                return
            endif
!            print *, "verticalCoordinate ", qt%verticalCoordinate
        endif

        if (l28(P_SURFS)) then
            if (qt%coherent) then
                call allocate_test(qt%surfs, qt%nosurfs, 1, "qt%surfs", moduleName)
            else
                call allocate_test(qt%surfs, qt%nosurfs, qt%noInstances, "qt%surfs", moduleName)
            endif
            call PVMIDLUnpack(qt%surfs, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking surfs" )
                return
            endif
!            print *, "surfs ", qt%surfs
        endif

        if (l28(P_INSTANCELEN)) then
            call PVMIDLUnpack(qt%instancelen, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instancelen" )
                return
            endif
!            print *, "instancelen ", qt%instancelen
        endif

        if (l28(P_VALUE)) then
            call Allocate_Test ( values, qt%instanceLen, qt%noInstances, &
            & 'values', ModuleName )
            call PVMIDLUnpack ( values, info, doreshape=.true. )
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking values" )
                return
            endif
!            print *, "values ", values
        endif

        ! need to fix this later to match the call to values
        if (l28(P_MASK)) then
            call CreateMaskArray ( mask, values )
            call PVMIDLUnpack(mask, info, doreshape=.true.)
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking mask" )
                return
            endif
!            print *, "mask ", mask
        endif

        ! make sure all hgrid-related fields are of the same shape, or at least allocated
        if (all(hshape /= (/0,0/))) then
            if (.not. l28(P_PHI)) then
                call Allocate_test( qt%phi, hshape(1), hshape(2), 'qt%phi', moduleName)
                qt%phi = 0.0_r8
            endif
            if (.not. l28(P_GEODLAT)) then
                call Allocate_test( qt%geodlat, hshape(1), hshape(2), 'qt%geodlat', moduleName)
                qt%geodlat = 0.0_r8
            endif
            if (.not. l28(P_LOSANGLE)) then
                call Allocate_test( qt%losangle, hshape(1), hshape(2), 'qt%losangle', moduleName)
                qt%losangle = 0.0_r8
            endif
            if (.not. l28(P_SOLARZENITH)) then
                call Allocate_test( qt%solarzenith, hshape(1), hshape(2), 'qt%solarzenith', moduleName)
                qt%solarzenith = 0.0_r8
            endif
            if (.not. l28(P_SOLARTIME)) then
                call Allocate_test( qt%solartime, hshape(1), hshape(2), 'qt%solartime', moduleName)
                qt%solartime = 0.0_r8
            endif
            if (.not. l28(P_TIME)) then
                call Allocate_test( qt%time, hshape(1), hshape(2), 'qt%time', moduleName)
                qt%time = 0.0_r8
            endif
            if (.not. l28(P_LONGITUDE)) then
                call Allocate_test( qt%lon, hshape(1), hshape(2), 'qt%lon', moduleName)
                qt%lon = 0.0_r8
            endif
        endif
    end subroutine

    subroutine ICFMSendVector (vec, tid, info)
        use QuantityPVM, only: PVMSendQuantity

        type(Vector_T), intent(in) :: vec
        integer, intent(in) :: tid
        integer, intent(out) :: info

        integer :: bufid
        integer :: numQty
        integer :: i

        call PVMFInitSend(PVMDataDefault, bufid)

        numQty = size(vec%quantities)
        call PVMIDLPack(numQty, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "packing numQty" )
            return
        endif

       do i = 1, numQty
            call PVMSendQuantity(vec%quantities(i), justpack=.true.)
        enddo

        call PVMFSend ( tid, QtyMsgTag, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector values" )
    end subroutine

    subroutine ICFMReceiveVector ( vec, qtydb, tid, callrecv)
        type(Vector_T), intent(out) :: vec
        type(QuantityTemplate_T), dimension(:), pointer :: qtydb
        integer, intent(in), optional :: TID ! Task to get it from
        logical, optional, intent(in) :: callrecv !true if this subroutine should call PVMFRecv, default is false

        integer :: info
        integer :: BUFFERID
        integer :: numQty
        integer :: i, j
        integer, dimension(:), pointer :: template

        ! Get buffer, we'll wait for it, assume the calling code knows it's coming.
        if (present(callrecv) .and. callrecv) call PVMFrecv ( tid, QtyMsgTag, bufferID )

        ! Now we unpack the information
        call PVMUnpackStringIndex ( vec%name, info )
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking name" )
            return
        endif
!        print *, "name ", vec%name

        call PVMIDLUnpack(numQty, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking number of quantities" )
            return
        endif
!        print *, "numQty ", numQty

        if (numQty == 0) return

        allocate(template(numQty), stat=info)
        if (info /= 0) call MLSMessage(MLSMSG_Error, ModuleName, "Cannot allocate template")

        allocate(vec%quantities(numQty), stat=info)
        if (info /= 0) call MLSMessage(MLSMSG_Error, ModuleName, "Cannot allocate vector quantities")

        do i=1, numQty
            call ICFMReceiveQuantity(vec%quantities(i)%template, vec%quantities(i)%values, vec%quantities(i)%mask)
            j = AddQuantityTemplateToDatabase(qtydb, vec%quantities(i)%template)
            template(i) = j
        end do

        vec%template = CreateVectorTemplate(qtydb, template)
        deallocate(template)

    end subroutine ICFMReceiveVector

end module

! $Log$
! Revision 1.1  2011/03/15 15:23:51  honghanh
! Initial imports
!
