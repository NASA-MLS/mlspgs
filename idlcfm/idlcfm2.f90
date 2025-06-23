module IDLCFM2_m
    use Allocate_Deallocate, only: Deallocate_test
    use CFM
    use CFM, only: QUANTITYTEMPLATE_T
    use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
    use PVM, only: PVMFRECV, PVMFINITSEND, PVMDATADEFAULT
    use MorePVM, only: PVMUnpackStringIndex, PVMUnpackLitIndex
    use MLSMessageModule, only: PVMERRORMESSAGE
    use QuantityTemplates, only: SetupNewQuantityTemplate, CopyQuantityTemplate
    use VectorsModule, only: RemapVectorMask, RemapVectorValue

    implicit none
    private
    public :: ICFMReceiveQuantity, ICFMReceiveVector, ICFMSendVector
!---------------------------- RCS Ident Info -------------------------------
    character(len=*), private, parameter :: ModuleName= &
        "$RCSfile$"
    private :: not_used_here
!---------------------------------------------------------------------------

    ! Local parameters
    integer, parameter, public :: QTYMSGTAG = 200

    integer, parameter, public :: SIG_SETUP   = 0
    integer, parameter, public :: SIG_CLEANUP = 1
    integer, parameter, public :: SIG_FWDMDL  = 2
    integer, parameter, public :: SIG_VECTOR  = 3
    integer, parameter, public :: SIG_DIE     = 4

    ! These methods are to help with debugging effort
    integer, parameter :: P_NAME                    = 1
    integer, parameter :: P_TYPE                    = 2
    integer, parameter :: P_OFFSET                  = 3
    integer, parameter :: P_COHERENT                = 4
    integer, parameter :: P_STACKED                 = 5
    integer, parameter :: P_LOGBASIS                = 6
    integer, parameter :: P_MINVALUE                = 7
    integer, parameter :: P_BADVALUE                = 8
    integer, parameter :: P_MOLECULE                = 9
    integer, parameter :: P_MODULE                  = 10
    integer, parameter :: P_SIGNAL                  = 11
    integer, parameter :: P_RADIOMETER              = 12
    integer, parameter :: P_FREQUENCIES             = 13
    integer, parameter :: P_FCOORD                  = 14
    integer, parameter :: P_NOCHANS                 = 15
    integer, parameter :: P_SURFS                   = 16
    integer, parameter :: P_VCOORD                  = 17
    integer, parameter :: P_NOSURFS                 = 18
    integer, parameter :: P_PHI                     = 19
    integer, parameter :: P_GEODLAT                 = 20
    integer, parameter :: P_LONGITUDE               = 21
    integer, parameter :: P_LOSANGLE                = 22
    integer, parameter :: P_SOLARZENITH             = 23
    integer, parameter :: P_SOLARTIME               = 24
    integer, parameter :: P_TIME                    = 25
    integer, parameter :: P_NOPROFS                 = 26
    integer, parameter :: P_INSTANCELEN             = 27
    integer, parameter :: P_VALUE                   = 28
    integer, parameter :: P_MASK                    = 29

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
        logical :: l29(29)
        type(QuantityTemplate_T) :: template

        if (present(callsend) .and. .not. callsend) mycallsend = .false.

        if (mycallsend .and. .not. present(tid)) then
            call output (MLSMSG_Error, "Missing 'tid'")
            return
        endif

        if (mycallsend) call PVMFInitSend(PVMDataDefault, bufid)

        template = qty%template

        l29(P_NAME) = template%name == 0
        l29(P_TYPE) = template%quantityType == 0
        l29(P_OFFSET) = .true.
        l29(P_COHERENT) = .true.
        l29(P_LOGBASIS) = .true.
        l29(P_MINVALUE) = template%logBasis !because if logBasis is false, no use for minValue
        l29(P_BADVALUE) = .true.
        l29(P_MOLECULE) = template%molecule == 0
        l29(P_MODULE) = template%instrumentModule == 0
        l29(P_SIGNAL) = template%signal == 0
        l29(P_RADIOMETER) = template%radiometer == 0
        l29(P_FREQUENCIES) = associated(template%frequencies) .and. size(template%frequencies) .gt. 0
        l29(P_FCOORD) = l29(P_FREQUENCIES)
        l29(P_NOCHANS) = template%noChans == 0
        l29(P_NOSURFS) = template%noSurfs == 0
        l29(P_SURFS) = allocated(template%surfs) .and. size(template%surfs) .gt. 0
        l29(P_VCOORD) = l29(P_SURFS)
        l29(P_PHI) = allocated(template%phi) .and. size(template%phi) .gt. 0
        l29(P_GEODLAT) = allocated(template%geodlat) .and. size(template%geodlat) .gt. 0
        l29(P_LONGITUDE) = allocated(template%lon) .and. size(template%lon) .gt. 0
        l29(P_LOSANGLE) = associated(template%losAngle) .and. size(template%losangle) .gt. 0
        l29(P_SOLARZENITH) = associated(template%solarZenith) .and. size(template%solarzenith) .gt. 0
        l29(P_SOLARTIME) = associated(template%solarTime) .and. size(template%solartime) .gt. 0
        l29(P_TIME) = associated(template%time) .and. size(template%time) .gt. 0
        ! this variable reflect the size of multiple arrays, so just send it for simplicity
        l29(P_NOPROFS) = .true.
        l29(P_INSTANCELEN) = .true.
        l29(P_VALUE) = associated(qty%values) .and. size(qty%values) .gt. 0
        l29(P_MASK) = associated(qty%mask) .and. size(qty%mask) .gt. 0

        call PVMIDLPack(l29, info)
        print *, "packing l29 ", l29

        if (mycallsend) then
            call PVMFSend ( tid, QtyMsgTag, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector value" )
        endif
    end subroutine ICFMSendQuantity

    ! subroutine ICFMReceiveQuantity ( QT, values, mask, tid, callrecv)
    subroutine ICFMReceiveQuantity ( QT, value1, mask1, tid, callrecv)
        use MLSCommon, only: R8
        use QuantityTemplates, only: QUANTITYTEMPLATE_T
        use ConstructQuantityTemplates, only: firstProperty, lastProperty, &
                                              propertyTable, p_majorFrame, &
                                              p_minorFrame, unitstable
        use parse_signal_m, only: parse_signal
        use MLSSignals_m, only: GetRadiometerIndex, GetModuleFromRadiometer, &
                                GetModuleIndex, GetRadiometerFromSignal, &
                                GetModuleFromSignal

        type (QuantityTemplate_T), intent(out) :: QT ! Template for quantity
        ! It's not inout, because then setupNewQuantityTemplate would deallocate
        ! the pointer components.  But the actual argument is put into a database
        ! using a shallow copy, so cleaning it up would clobber a database item.
        real (r8), dimension(:), pointer :: VALUE1 ! Values for quantity
        character, dimension(:), pointer :: MASK1 ! Mask
        integer, intent(in), optional :: TID ! Task to get it from
        logical, optional, intent(in) :: callrecv !true if this subroutine should call PVMFRecv, default is false

        ! Local variables
        integer :: BUFFERID                 ! From pvm
        integer :: INFO                     ! Flag
        integer, parameter :: P_LAST = P_MASK + 1
        logical :: l29(P_LAST - 1)
        logical :: PROPERTIES(firstProperty : lastProperty) ! Properties for this quantity type
        character(len=640) :: signalString
        integer, dimension(:), pointer :: SignalInds ! From parse signal
        logical, pointer :: channels(:)     ! From Parse_Signal
        integer :: sideband
        character(len=16) :: radiometerString
        integer, dimension(2) :: hshape ! shape of hgrid-related fields

        ! Executable code
        ! First, sanitize our input
        if (associated(value1) .or. associated(mask1)) then
            call output(MLSMSG_Error, &
            "'values' and 'mask' parameters must be nullified before receiving quantity")
            return
        endif

        hshape = 0
        call InitializeQuantityTemplate(qt)
        nullify(channels, signalinds)

        ! Get buffer, we'll wait for it, assume the calling code knows it's coming.
        if (present(callrecv)) then
            if(callrecv) call PVMFrecv ( tid, QtyMsgTag, bufferID )
        endif

        call PVMIDLUnpack(l29, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking l29." )
            call clearout
            return
        endif
!        print *, "l29 ", l29

        if (l29(P_NAME)) then
            call PVMUnpackStringIndex ( qt%name, info)
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking name" )
                call clearout
                return
            endif
!            print *, "name " , qt%name
        endif

        if (l29(P_TYPE)) then
            call PVMUnpackLitIndex ( qt%quantityType, info )
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking quantityType" )
                call clearout
                return
            endif
!            print *, "type " , qt%quantityType
            if (qt%quantityType < 0) then
                call output(MLSMSG_Error, "invalid quantity type")
                call clearout
                return
            endif
            qt%unit = unitsTable(qt%quantityType)
!            print *, "unit ", qt%unit
            properties = propertyTable(:, qt%quantityType)
            qt%majorFrame = properties(p_majorFrame)
!            print *, "majorFrame ", qt%majorFrame
            qt%minorFrame = properties(p_minorFrame)
!            print *, "minorFrame ", qt%minorFrame
        endif

        if (l29(P_OFFSET)) then
            call PVMIDLUnpack(qt%instanceOffset, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instanceOffset" )
                call clearout
                return
            endif
!            print *, "instanceOffset ", qt%instanceOffset
        endif

        if (l29(P_COHERENT)) then
            call PVMIDLUnpack(qt%coherent, info)
            if (info /= 0) then
                call PVMErrorMessage(info, "unpacking coherent")
                call clearout
                return
            endif
!            print *, "coherence ", qt%coherent
        endif

        if (l29(P_STACKED)) then
            call PVMIDLUnpack(qt%stacked, info)
            if (info /= 0) then
                call PVMErrorMessage(info, "unpacking stacked")
                call clearout
                return
            endif
!            print *, "stack ", qt%stack
        endif

        if (l29(P_LOGBASIS)) then
            call PVMIDLUnpack(qt%logbasis, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking logbasis" )
                call clearout
                return
            endif
!            print *, "logbasis ", qt%logbasis
        endif

        if (l29(P_MINVALUE)) then
            call PVMIDLUnpack(qt%minvalue, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking minvalue" )
                call clearout
                return
            endif
!            print *, "minvalue ", qt%minvalue
        endif

        if (l29(P_BADVALUE)) then
            call PVMIDLUnpack(qt%badvalue, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking badvalue" )
                call clearout
                return
            endif
!            print *, "badvalue ", qt%badvalue
        endif

        if (l29(P_MOLECULE)) then
            call PVMUnpackLitIndex(qt%molecule, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking molecule" )
                call clearout
                return
            endif
!            print *, "molecule ", qt%molecule
        endif

        if (l29(P_MODULE)) then
            call PVMIDLUnpack(signalString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instrumentModule" )
                call clearout
                return
            endif
            call GetModuleIndex(signalString, qt%instrumentModule)
!            print *, "instrumentModule ", qt%instrumentModule
        endif

        if (l29(P_SIGNAL)) then
            call PVMIDLUnpack(signalString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking signalString" )
                call clearout
                return
            endif
!            print *, "signalString ", signalString
            call parse_Signal ( signalString, signalInds, sideband=sideband, channels=channels)
            qt%signal = signalInds(1)
            qt%sideband = sideband
            call deallocate_test ( signalInds, 'signalInds', ModuleName )
!            print *, "signal ", qt%signal
            qt%radiometer = GetRadiometerFromSignal(qt%signal)
            qt%instrumentModule = GetModuleFromSignal(qt%signal)
        endif

        if (l29(P_RADIOMETER)) then
            call PVMIDLUnpack(radiometerString, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking radiometerString" )
                call clearout
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

        if (l29(P_NOCHANS)) then
            call PVMIDLUnpack(qt%noChans, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noChans" )
                call clearout
                return
            endif
!            print *, "noChans ", qt%noChans
        endif

        if (l29(P_FCOORD)) then
            call PVMUnpackLitIndex(qt%frequencyCoordinate, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking frequencyCoordinate" )
                call clearout
                return
            endif
!            print *, "frequencyCoordinate ", qt%frequencyCoordinate
        endif

        if (l29(P_FREQUENCIES)) then
            allocate(qt%frequencies(qt%noChans), stat=info)
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
            endif
            call PVMIDLUnpack(qt%frequencies, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking frequencies" )
                call clearout
                return
            endif
!            print *, "frequencies ", qt%frequencies
        endif

        if (l29(P_NOSURFS)) then
            call PVMIDLUnpack(qt%noSurfs, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noSurfs" )
                call clearout
                return
            endif
!            print *, "noSurfs ", qt%noSurfs
        endif

        if (l29(P_NOPROFS)) then
            call PVMIDLUnpack(qt%noInstances, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking noInstances" )
                call clearout
                return
            endif
!            print *, "noProfs ", qt%noInstances
        endif

        if (l29(P_PHI)) then
            if (qt%stacked) then
                allocate(qt%phi(1, qt%noInstances), stat=info)
            else
                allocate(qt%phi(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output (MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%phi, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking phi" )
                call clearout
                return
            endif
!            print *, "phi ", qt%phi
            hshape = shape(qt%phi)
            ! None of the idl procedures know about cross angles yet, so they
            ! won't be sending us any news about that geolocation until
            ! someone goes in and educates the idl procedures
            ! in 
            ! mlspgs/idlcfm/idl
            ! For now we resort to the following hackery-quackery
            allocate(qt%crossAngles(1), stat=info)
            if (info /= 0) then
                call output (MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            qt%crossAngles = 0.
            qt%noCrossTrack = 1
        endif

        if (l29(P_GEODLAT)) then
            if (qt%stacked) then
                allocate(qt%geodlat(1, qt%noInstances), stat=info)
            else
                allocate(qt%geodlat(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output (MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%geodlat, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking geodlat" )
                call clearout
                return
            endif
!            print *, "geodlat ", qt%geodlat
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%geodlat))) then
                call output(MLSMSG_Warning, "geodlat's shape is not the same as others' shape")
            endif
            hshape = shape(qt%geodlat)
        endif

        if (l29(P_LONGITUDE)) then
            if (qt%stacked) then
                allocate(qt%lon(1, qt%noInstances), stat=info)
            else
                allocate(qt%lon(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%lon, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking lon" )
                call clearout
                return
            endif
!            print *, "lon ", qt%lon
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%lon))) then
                call output(MLSMSG_Warning, "lon's shape is not the same as others' shape")
            endif
            hshape = shape(qt%lon)
        endif

        if (l29(P_LOSANGLE)) then
            if (qt%stacked) then
                allocate(qt%losAngle(1, qt%noInstances), stat=info)
            else
                allocate(qt%losAngle(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output (MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%losAngle, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking losAngle" )
                call clearout
                return
            endif
!            print *, "losAngle ", qt%losAngle
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%losAngle))) then
                call output(MLSMSG_Warning, "losAngle's shape is not the same as others' shape")
            endif
            hshape = shape(qt%losAngle)
        endif

        if (l29(P_SOLARZENITH)) then
            if (qt%stacked) then
                allocate(qt%solarZenith(1, qt%noInstances), stat=info)
            else
                allocate(qt%solarZenith(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%solarZenith, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking solarZenith" )
                call clearout
                return
            endif
!            print *, "solarZenith ", qt%solarZenith
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%solarZenith))) then
                call output(MLSMSG_Warning, "solarZenith's shape is not the same as others' shape")
            endif
            hshape = shape(qt%solarZenith)
        endif

        if (l29(P_SOLARTIME)) then
            if (qt%stacked) then
                allocate(qt%solartime(1, qt%noInstances), stat=info)
            else
                allocate(qt%solartime(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%solartime, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking solartime" )
                call clearout
                return
            endif
!            print *, "solartime ", qt%solartime
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%solartime))) then
                call output(MLSMSG_Warning, "solarTime's shape is not the same as others' shape")
            endif
            hshape = shape(qt%solartime)
        endif

        if (l29(P_TIME)) then
            if (qt%stacked) then
                allocate(qt%time(1, qt%noInstances), stat=info)
            else
                allocate(qt%time(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%time, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking time" )
                call clearout
                return
            endif
!            print *, "time ", qt%time
            if (all(hshape /= (/0,0/)) .and. all(hshape /= shape(qt%time))) then
                call output(MLSMSG_Warning, "time's shape is not the same as others' shape")
            endif
            hshape = shape(qt%time)
        endif

        if (l29(P_VCOORD)) then
            call PVMUnpackLitIndex(qt%verticalCoordinate, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking verticalCoordinate" )
                call clearout
                return
            endif
!            print *, "verticalCoordinate ", qt%verticalCoordinate
        endif

        if (l29(P_SURFS)) then
            if (qt%coherent) then
                allocate(qt%surfs(qt%nosurfs, 1), stat=info)
            else
                allocate(qt%surfs(qt%nosurfs, qt%noInstances), stat=info)
            endif
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            call PVMIDLUnpack(qt%surfs, info, doreshape=.true.)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking surfs" )
                call clearout
                return
            endif
!            print *, "surfs ", qt%surfs
        endif

        if (l29(P_INSTANCELEN)) then
            call PVMIDLUnpack(qt%instancelen, info)
            if (info /=0 ) then
                call PVMErrorMessage ( info, "unpacking instancelen" )
                call clearout
                return
            endif
!            print *, "instancelen ", qt%instancelen
            if (qt%instancelen == 0) then
                call output (MLSMSG_Warning, "Quantity template has instanceLen 0")
            endif 
        endif

        if (l29(P_VALUE)) then
            allocate( value1(qt%instanceLen*qt%noInstances), stat=info)
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of machine")
                call clearout
                return
            endif
            call PVMIDLUnpack ( value1, info )
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking values" )
                call clearout
                return
            endif
!            print *, "values ", values
        endif

        ! need to fix this later to match the call to values
        if (l29(P_MASK)) then
            allocate(mask1(size(value1)), stat=info)
            if (info /= 0) then
                call output(MLSMSG_Error, "Out of memory")
                call clearout
                return
            endif
            mask1 = char(0) ! All vector elements are interesting
            call PVMIDLUnpack(mask1, info)
            if ( info /= 0 ) then
                call PVMErrorMessage ( info, "unpacking mask" )
                call clearout
                return
            endif
!            print *, "mask ", mask
        endif

        ! make sure all hgrid-related fields are of the same shape, or at least allocated
        if (all(hshape /= (/0,0/))) then
            if (.not. l29(P_PHI)) then
                allocate( qt%phi(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%phi = 0.0_r8
            endif
            if (.not. l29(P_GEODLAT)) then
                allocate( qt%geodlat(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%geodlat = 0.0_r8
            endif
            if (.not. l29(P_LOSANGLE)) then
                allocate( qt%losangle(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%losangle = 0.0_r8
            endif
            if (.not. l29(P_SOLARZENITH)) then
                allocate( qt%solarzenith(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%solarzenith = 0.0_r8
            endif
            if (.not. l29(P_SOLARTIME)) then
                allocate( qt%solartime(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%solartime = 0.0_r8
            endif
            if (.not. l29(P_TIME)) then
                allocate( qt%time(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%time = 0.0_r8
            endif
            if (.not. l29(P_LONGITUDE)) then
                allocate( qt%lon(hshape(1), hshape(2)), stat=info)
                if (info /= 0) then
                    call output(MLSMSG_Error, "Out of memory")
                    call clearout
                    return
                endif
                qt%lon = 0.0_r8
            endif
        endif

        contains
        
        subroutine clearout
            if(associated(channels)) deallocate(channels) 

            if(associated(qt%frequencies)) then
                deallocate(qt%frequencies)
                nullify(qt%frequencies)
            endif

            if (allocated(qt%phi)) then
                deallocate(qt%phi)
            endif

            if (associated(qt%crossAngles)) then
                deallocate(qt%crossAngles)
            endif

            if (allocated(qt%geodlat)) then
                deallocate(qt%geodlat)
            endif

            if (allocated(qt%lon)) then
                deallocate(qt%lon)
            endif

            if (associated(qt%losangle)) then
                deallocate(qt%losangle)
                nullify(qt%losangle)
            endif

            if (associated(qt%solarZenith)) then
                deallocate(qt%solarzenith)
                nullify(qt%solarzenith)
            endif

            if (associated(qt%solartime)) then
                deallocate(qt%solartime)
                nullify(qt%solartime)
            endif

            if (associated(qt%time)) then
                deallocate(qt%time)
                nullify(qt%time)
            endif

            if (allocated(qt%surfs)) then
                deallocate(qt%surfs)
            endif

            if (associated(value1)) then
                deallocate(value1)
                nullify(value1)
            endif

            if (associated(mask1)) then
                deallocate(mask1)
                nullify(mask1)
            endif
        end subroutine clearout
    end subroutine ICFMReceiveQuantity

    subroutine InitializeQuantityTemplate (template)
        ! This function is to initialize a template
        ! to what IDLCFM library need
        type(QuantityTemplate_T), intent(out) :: template

        template%name = 0
        template%quantityType = 0
        template%noInstances = 1
        template%noSurfs = 1
        template%noChans = 1
        template%NoCrossTrack = 1
        template%coherent = .true.
        template%stacked = .true.
        template%regular = .true. ! we don't have irregular quantities
        template%minorFrame = .false.
        template%majorFrame = .false.
        template%logBasis = .false.
        template%minValue = -huge(0.0_r8)
        template%noInstancesLowerOverlap = 0
        template%noInstancesUpperOverlap = 0
        template%badValue = huge(0.0_r8)
        template%unit = 0
        template%instanceLen = 0
        template%verticalCoordinate = l_none
        template%sharedVGrid = .false.
        template%vGridIndex = 0
        nullify(template%time, template%solartime, template%crossAngles)
        nullify(template%solarzenith, template%losangle)
        nullify(template%chaninds, template%channels, template%frequencies)
        template%fgridindex = 0
        template%frequencyCoordinate = l_none
        template%lo = 0.0_r8
        template%sharedfgrid = .false.
        template%sideband = 0
        template%signal = 0
        template%instrumentModule = 0
        template%radiometer = 0
        template%reflector = 0
        template%molecule = 0
        nullify(template%surfindex, template%chanIndex)
    end subroutine InitializeQuantityTemplate

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
    end subroutine ICFMSendVector

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
        integer :: type

        ! Get buffer, we'll wait for it, assume the calling code knows it's coming.
        if (present(callrecv)) then
            if (callrecv) then
                call PVMFrecv ( tid, QtyMsgTag, bufferID )
                call PVMIDLUnpack(type, info)
                if (info /= 0) then
                    call PVMErrorMessage(info, "unpacking type")
                    return
                endif

                if (type /= SIG_VECTOR) then
                    call output(MLSMSG_Error, "ICFMReceiveVector: the received is not vector")
                    return
                endif
            endif
        endif

        ! Now we unpack the information
        call PVMUnpackStringIndex ( vec%name, info )
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking name" )
            return
        endif
        print *, "name ", vec%name

        call PVMIDLUnpack(numQty, info)
        if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking number of quantities" )
            return
        endif
        print *, "numQty ", numQty

        if (numQty == 0) return

        allocate(template(numQty), stat=info)
        if (info /= 0) then
            call output (MLSMSG_Error, "Cannot allocate template")
            return
        endif

        allocate(vec%quantities(numQty), stat=info)
        if (info /= 0) then
            call output(MLSMSG_Error, "Cannot allocate vector quantities")
            deallocate(template)
            return
        endif

        do i=1, numQty
            call ICFMReceiveQuantity ( vec%quantities(i)%template, &
              & vec%quantities(i)%value1, vec%quantities(i)%mask1 )
            if ( associated(vec%quantities(i)%value1) ) &
              & call remapVectorValue ( vec%quantities(i) )
            if ( associated(vec%quantities(i)%mask1) ) &
              & call remapVectorMask ( vec%quantities(i) )
            vec%quantities(i)%index = i
            j = AddQuantityTemplateToDatabase(qtydb, vec%quantities(i)%template)
            template(i) = j
        end do

        vec%template = CreateVectorTemplate(qtydb, template)
        deallocate(template)

    end subroutine ICFMReceiveVector

    subroutine output (level, msg)
        integer, intent(in) :: level
        character(len=*), intent(in) :: msg

        select case (level)
        case (MLSMSG_Success)
            print *, "succeed: ", msg
        case (MLSMSG_Info)
            print *, "info: ", msg
        case (MLSMSG_Debug)
            print *, "debug: ", msg
        case (MLSMSG_Warning) 
            print *, "warning: ", msg
        case (MLSMSG_Error) 
            print *, "error: ", msg
        case (MLSMSG_Crash)
            print *, "crash: ", msg
        case default
            print *, "unknown: ", msg
        end select
    end subroutine output

end module

! $Log$
! Revision 1.11  2016/02/04 22:00:57  pwagner
! added SIG_DIE; define P_stuff in just one place now
!
! Revision 1.10  2016/01/07 17:56:06  pwagner
! Reflects new crossAngles, pointer remapping
!
! Revision 1.9  2015/08/17 18:21:27  pwagner
! Changed to reflect loss of pointiness by some qty template components
!
! Revision 1.8  2012/01/09 22:36:54  pwagner
! Workaround for ifc 12 bug
!
! Revision 1.7  2012/01/03 17:21:25  honghanh
! Remove unused variables, and incorporate changes from
! lit_parm and Molecules_M
!
! Revision 1.6  2011/09/07 06:34:46  honghanh
! Make modification to send matrix,and add more error handling
!
! Revision 1.5  2011/06/27 21:28:52  honghanh
! Fixed bug in a few if statement since Fortran does not have
! lazy evaluation of boolean expressions
!
! Revision 1.4  2011/05/27 06:12:48  honghanh
! Add warning message upon unpacking instanceLen 0
!
! Revision 1.3  2011/05/27 06:08:41  honghanh
! Add VectorHandler and type code of ICFMReceiveVector to receive vector independently of the call to ForwardModel
!
! Revision 1.2  2011/04/16 22:03:45  honghanh
! *** empty log message ***
!
! Revision 1.1  2011/03/15 15:23:51  honghanh
! Initial imports
!
