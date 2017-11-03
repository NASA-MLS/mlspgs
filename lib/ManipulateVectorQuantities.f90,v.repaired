! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ManipulateVectorQuantities ! Various routines for manipulating vectors

  ! This modules contains routines needed for manipulating vectors.

   use HighOutput, only: OutputnamedValue
   use MLSMessagemodule, only: MLSMSG_Error, MLSMessage
   use MLSKinds, only: R8, Rv
   use MLSNumerics, only: Hunt
   use Molecules, only: L_Rhi
   use Trace_M, only: Trace_Begin, Trace_End
   use Vectorsmodule, only: VectorValue_T, Vector_T, Dump
   use Intrinsic, only: L_Calsidebandfraction, L_Channel, L_Columnabundance, &
     & L_Isotoperatio, L_Limbsidebandfraction, L_None, L_Phitan, L_Radiance, &
     & L_Tscat, L_Vmr, Phyq_ProFiles

  implicit none

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  private

  public :: AnyGoodDataInQty, DoHGridsMatch, DoQuantitiesMatch, &
    & DoVGridsMatch, DoVGridsMatch_VEC, &
    & DoFGridsMatch, DoQtysDescribeSameThing, DoVectorsMatch, &
    & FillWithCombinedChannels, FindClosestInstances, FindOneClosestInstance, &
    & FindInstanceWindow

  interface DoVGridsMatch
    module procedure DoVGridsMatch_Vec
  end interface

  logical, parameter               :: DEEBUG = .false.

contains

  ! --------------------------------------- AnyGoodDataInQty --------------
  logical function AnyGoodDataInQty ( a, &
    & precision, quality, status, quality_min, op )
    ! Returns true if any of the mask != char(0)
    ! Returns true if any of the precision values >= 0 (if present)
    ! Returns true if any of the quality values >= quality_min (if present)
    ! Returns true if any of the status values are not odd (if present)
    !
    ! By default the result returned is the "and" of all the results
    ! for each of its args
    ! If instead you want the "or", set op = "or"
    
    ! If mask is not associated, should we return TRUE or FALSE?
    ! The answer should be TRUE
    ! But if values aren't associated for precision, quality, or status,
    ! the behavior depends on context and may be the opposite 
    ! of what you would expect
    ! If you supply only precision, and that precision array is not associated
    ! we return TRUE
    ! 
    ! Should we warn of misuse if user supplies quality w/o quality_min?
    use MLSStrings, only: LowerCase

    type ( VectorValue_T ), pointer   , optional :: a         ! qty values
    type ( VectorValue_T ), pointer   , optional :: precision ! Precision
    type ( VectorValue_T ), pointer   , optional :: quality   ! Quality
    type ( VectorValue_T ), pointer   , optional :: status    ! Status
    real(rv),               intent(in), optional :: quality_min
    character(len=*),       intent(in), optional :: op
    ! Local variables
    character(len=4) :: myOp

    ! Executable code
    myOp = "and"
    if ( present(op) ) myOp = lowerCase(op)
    AnyGoodDataInQty = .true.
    if ( myOp /= "and" ) AnyGoodDataInQty = .false.

    if ( present(a) ) then
      if ( associated(a) ) then
        AnyGoodDataInQty = .true. ! Return TRUE if mask not associated
        if ( associated ( a%mask ) ) &
        & call joinResults( any(a%mask == char(0) ) )
      endif
    endif
    if ( present(precision) ) then
      if ( associated(precision) ) then
        if ( associated ( precision%values ) ) &
        & call joinResults( any(precision%values >= 0._rv) )
      endif
    endif
    if ( present(quality) .and. present(quality_min) ) then
      if ( associated(quality) ) then
        if ( associated ( quality%values ) ) &
        & call joinResults( any(quality%values >= quality_min) )
      endif
    endif
    if ( present(status) ) then
      if ( associated(status) ) then
        if ( associated ( status%values ) ) &
        & call joinResults( any(mod(int(status%values), 2) /= 0) )
      endif
    endif
    
    contains
    subroutine joinResults (arg)
      logical, intent(in) :: arg
      if ( myOp == "and" ) then
        AnyGoodDataInQty = ( AnyGoodDataInQty .and. arg )
      else
        AnyGoodDataInQty = ( AnyGoodDataInQty .or. arg )
      endif
    end subroutine joinResults
    
  end function AnyGoodDataInQty

  ! ------------------------------ FindClosestInstances -----------------
  subroutine FindClosestInstances ( referenceQuantity, soughtQuantity,&
    referenceIndices )
    ! This subroutine is similar to FindOneClosestInstance (and calls it in
    ! fact), and finds an array of closest instances

    ! Dummy arguments
    type (VectorValue_T), intent(in) :: referenceQuantity ! e.g. temperature
    type (VectorValue_T), intent(in) :: soughtQuantity ! e.g. ptan, radiance
    integer, dimension(soughtQuantity%template%noInstances), &
      & intent(out) :: referenceIndices

    ! Local variables
    integer :: instance

    ! Executable code
    do instance = 1, soughtQuantity%template%noInstances
      referenceIndices(instance) = FindOneClosestInstance ( &
        & referenceQuantity, soughtQuantity, instance )
    end do
  end subroutine FindClosestInstances

  ! ---------------------------------------- FindOneClosestInstance -----
  integer function FindOneClosestInstance ( referenceQuantity, &
    soughtQuantity, instance )
    use hGridsDatabase, only: FindClosestMatch
    use Intrinsic, only: L_Time
    ! This returns the instance index into a stacked quantity for the
    ! instance 'closest' to the given instance in an unstacked one
    type (VectorValue_T), intent(in) :: ReferenceQuantity ! e.g. temperature
    type (VectorValue_T), intent(in), target :: SoughtQuantity ! e.g. ptan, radiance
    integer, intent(in) :: Instance

    ! Local variables
    integer :: horizontalCoordinate
    real (r8), dimension(:,:), pointer :: Seek ! The thing to look for

    ! Executable:
    horizontalCoordinate = referenceQuantity%template%horizontalCoordinate
    ! We'll skip the error checking we could do at this point, for speed.

    ! First we'll do a hunt to get ourselves in the right area.  Might as
    ! well start looking where we think that will be.

    if ( soughtQuantity%template%quantityType == l_phiTan ) then
      seek => soughtQuantity%values
    else if ( horizontalCoordinate == l_time ) then
      seek => soughtQuantity%template%time
    else
      seek => soughtQuantity%template%phi
    end if
    ! Call FindClosestMatch to do the work
    if ( horizontalCoordinate == l_time ) then
      FindOneClosestInstance = FindClosestMatch ( &
        & referenceQuantity%template%time(1,:), &
        & seek, instance )
    else if ( allocated(referenceQuantity%template%phi) ) then
      FindOneClosestInstance = FindClosestMatch ( &
        & referenceQuantity%template%phi(1,:), &
        & seek, instance )
    else
      call Dump ( referenceQuantity%template )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'In FindOneClosestInstance, reference phi is not allocated and ' // &
        & 'horizontal coordinate is not time' )
    end if
  end function FindOneClosestInstance

  ! --------------------------------------- FindInstanceWindow ---------
  subroutine FindInstanceWindow ( Quantity, PhiTan, MAF, PhiWindow, &
    & WindowUnits, WindowStart, WindowFinish )
    ! This returns the start and end of a window into a quantity such as
    ! temperature for a given instance of a minor frame quantity
    type (VectorValue_T), intent(in) :: Quantity ! Quantity e.g. temperature
    type (VectorValue_T), intent(in) :: PhiTan ! Phitan information
    integer, intent(in) :: MAF            ! Major frame sought
    real (r8), intent(in) :: PhiWindow(2) ! Window size before and after PhiTan
    integer, intent(in) :: WindowUnits    ! PHYQ_Angle or PHYQ_Profiles
    integer, intent(out) :: WindowStart   ! Output window start
    integer, intent(out) :: WindowFinish  ! Output window finish

    ! Internal variables
    integer :: CLOSESTINSTANCE
    integer :: Me = -1                  ! String index for trace cacheing
    real(r8) :: PHIMIN, PHIMAX          ! Limiting values of phi for this MAF

    ! Executable code
    call trace_begin ( me, 'FindInstanceWindow', cond=.false. )
    ! WindowUnits was checked to be either PHYQ_Profiles or PHYQ_Angle.
    ! If WindowUnits had been given as 0 or 0:0, the closest instance will be
    ! chosen.
    if ( windowUnits == PHYQ_Profiles ) then
      ! Set the window start : window finish so that there are phiWindow(1) 
      ! profiles before the closest instance, and phiWindow(2) after the
      ! closest instance.
      closestInstance = FindOneClosestInstance ( quantity, phiTan, maf )
      windowStart = min ( quantity%template%noInstances, &
        & max ( 1, closestInstance - nint ( phiWindow(1) ) ) )
      windowFinish = min ( quantity%template%noInstances, &
        & closestInstance + nint ( phiWindow(2) ) )
    else ! windowUnits == PHYQ_Angle
      phiMin = minval ( phiTan%values(:,maf) ) - phiWindow(1)
      phiMax = maxval ( phiTan%values(:,maf) ) + phiWindow(2)
      call Hunt ( quantity%template%phi(1,:), phiMin, windowStart, &
        & allowTopValue=.true. )
      call Hunt ( quantity%template%phi(1,:), phiMax, windowFinish, &
        & allowTopValue=.true. )
      windowStart = min ( quantity%template%noInstances, &
        & max ( 1, windowStart - 1 ) )
      windowFinish = min ( quantity%template%noInstances, &
        & max ( 1, windowFinish + 1 ) )
    end if
    call trace_end ( cond=.false. )
  end subroutine FindInstanceWindow

  ! -------------------------------------- FillWithCombinedChannels ----------
  subroutine FillWithCombinedChannels ( quantity, sourceQuantity, message, mapping )
    use Matrixmodule_0, only: MatrixElement_T, M_Full, Createblock, Sparsify
    use MLSSignals_M, only: Signal_T, Getsignal
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    ! This routine takes a (typically radiance) quantity on one set of channels
    ! and combines the channels together appropriately to make them representative
    ! of the data in another (presumably at a finer resolution.
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Quantitiy on finer(?) channel grid
    character (len=*), intent(out) :: MESSAGE ! Possible error message
    type (MatrixElement_T), intent(inout), optional :: MAPPING ! A matrix_0 mapping
    ! Local variables
    type (Signal_T) :: signal, sourceSignal
    integer :: COUT, CIN              ! Channel counters
    integer :: IOUT, IIN              ! Indices
    integer :: SURF                   ! Loop counter
    real(r8) :: CENTER, HALFWIDTH     ! Channel locations
    integer, dimension(:), pointer :: NOINSIDE ! Number of channels that were caught

    ! Do some sanity checking
    message = ''
    if ( .not. DoVGridsMatch ( quantity, sourceQuantity ) ) then
      message = 'Quantities must have matching vGrids'
      return
    end if
    if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) ) then
      message = 'Quantities must have matching hGrids'
      return
    end if

    if ( quantity%template%signal == 0 .or. &
      &  quantity%template%frequencyCoordinate /= l_channel ) then
      message = 'Quantity must have channels'
      return
    end if
    if ( sourceQuantity%template%signal == 0 .or. &
      &  sourceQuantity%template%frequencyCoordinate /= l_channel ) then
      message = 'source quantity must have channels'
      return
    end if

    signal = GetSignal ( quantity%template%signal )
    sourceSignal = GetSignal ( sourceQuantity%template%signal )

    if ( signal%radiometer /= sourceSignal%radiometer ) then
      message = 'quantities must be from same radiometer'
      return
    end if

    nullify ( noInside )
    call Allocate_test ( noInside, quantity%template%noChans, 'noInside', ModuleName )

    ! Possibly setup a mapping matrix
    if ( present ( mapping ) ) then
      call CreateBlock ( mapping, &
        & quantity%template%instanceLen, sourceQuantity%template%instanceLen, &
        & kind=m_full, nChan=quantity%template%noChans )
    end if

    ! Do a quick first pass to get the numbers in each output channel
    do cOut = 1, quantity%template%noChans
      ! Work out where this output channel is
      center = signal%centerFrequency + signal%direction * signal%frequencies ( cOut )
      halfWidth = signal%widths ( cOut ) / 2.0
      ! Now map this back into sourceSignal%frequencies
      center = ( center - sourceSignal%centerFrequency ) * sourceSignal%direction
      noInside(cOut) = count ( &
        & ( sourceSignal%frequencies >= ( center - halfWidth ) ) .and. &
        & ( sourceSignal%frequencies < ( center + halfWidth ) ) )
    end do
    if ( any ( noInside == 0 ) ) then
      message = 'Some channels have no source'
      return
    end if

    ! Loop over the channels in the result
    do cOut = 1, quantity%template%noChans
      ! Work out where this output channel is
      center = signal%centerFrequency + signal%direction * signal%frequencies ( cOut )
      halfWidth = signal%widths ( cOut ) / 2.0
      ! Now map this back into sourceSignal%frequencies
      center = ( center - sourceSignal%centerFrequency ) * sourceSignal%direction
      ! Use nested loops to make the code easier to read, rather than complicated indexing
      do cIn = 1, sourceQuantity%template%noChans
        if ( ( sourceSignal%frequencies(cIn) >= ( center - halfWidth ) ) .and. &
          &  (  sourceSignal%frequencies(cIn) < ( center + halfWidth ) ) ) then
          do surf = 0, quantity%template%noSurfs - 1 ! 0..n-1 makes indexing easier
            iOut = cOut + surf * quantity%template%noChans
            iIn = cIn + surf * sourceQuantity%template%noChans
            quantity%values ( iOut, : ) = quantity%values ( iOut, : ) + &
              & sourceQuantity%values ( iIn, : ) / noInside(cOut)
            ! Possibly fill mapping matrix
            if ( present ( mapping ) ) then
              mapping%values ( iOut, iIn ) = 1.0 / noInside(cOut)
            end if
          end do
        end if
      end do
    end do
    call Sparsify ( mapping )

    call Deallocate_test ( noInside, 'noInside', ModuleName )
  end subroutine FillWithCombinedChannels

  ! --------------------------------------- DoHGridsMatch --------------
  logical function DoHGridsMatch ( a, b, spacingOnly )
    ! Returns true if quantities have same hGrid information
    type (VectorValue_T), intent(in) :: A ! First quantity
    type (VectorValue_T), intent(in) :: B ! Second quantity
    logical, optional, intent(in) :: SPACINGONLY

    ! Local parameters
    logical :: MYSPACINGONLY
    real (r8), parameter :: PHITOL = 0.01 ! Tolerance in angle
    real (r8) :: MINA, MINB, MAXA, MAXB ! Information on a and b

    ! Executable code
    DoHGridsMatch = .false.

    if ( associated(a%template%the_hGrid) .and. &
       & associated(b%template%the_hGrid) ) then ! How could this fail?
      if ( a%template%the_hGrid%type /= b%template%the_hGrid%type ) return
    end if

    mySpacingOnly = .false.
    if ( present ( spacingOnly ) ) mySpacingOnly = spacingOnly

    if ( a%template%stacked .neqv. b%template%stacked ) return
    if ( mySpacingOnly .and. .not. a%template%stacked ) return

    if ( .not. mySpacingOnly ) then
      if ( a%template%noInstances /= b%template%noInstances ) return
      if ( .not. ( allocated(a%template%phi) .eqv. &
                 & allocated(b%template%phi) ) ) return
      if ( allocated(a%template%phi) ) then
        if ( any(abs(a%template%phi - &
          &          b%template%phi) > PhiTol) ) return
      end if
      DoHGridsMatch = .true.
    else
      ! Here we default to true
      doHGridsMatch = .true.
      if ( a%template%noInstances == 1 ) return
      if ( b%template%noInstances == 1 ) return
      mina = minval ( &
        & a%template%phi(1,2:a%template%noInstances) - &
        & a%template%phi(1,1:a%template%noInstances-1) )
      minb = minval ( &
        & b%template%phi(1,2:b%template%noInstances) - &
        & b%template%phi(1,1:b%template%noInstances-1) )
      maxa = maxval ( &
        & a%template%phi(1,2:a%template%noInstances) - &
        & a%template%phi(1,1:a%template%noInstances-1) )
      maxb = minval ( &
        & b%template%phi(1,2:b%template%noInstances) - &
        & b%template%phi(1,1:b%template%noInstances-1) )
      doHGridsMatch = all ( (/ &
        & maxa-mina, maxb-minb, abs(maxa-maxb), abs(mina-minb) /) < phiTol )
    end if

  end function DoHGridsMatch

  ! ------------------------------------- DoQuantitiesMatch --
  logical function DoQuantitiesMatch ( a, b, options )
    ! Returns true if quantities share all important template information.
    type (VectorValue_T), intent(in) :: A ! First quantity
    type (VectorValue_T), intent(in) :: B ! Second quantity
    character(len=*), optional, intent(in)      :: options ! e.g. 's'[trict]
    ! logical, optional, intent(in)      :: Strict ! Must every attribute match?
    !                                   Or just the ones we actually write out
    logical :: myStrict
    integer, parameter :: NSCREENEDTYPES = 8
    integer, dimension(NSCREENEDTYPES) :: SCREENEDTYPES
    ! Executable code
    screenedtypes = (/ &
      & l_vmr, l_columnabundance, l_rhi, l_isotoperatio, l_radiance, l_calsidebandfraction, l_limbsidebandfraction, l_Tscat &
      & /)
    myStrict = .false.
    if ( present(options) ) myStrict = index(options, 's') > 0
    DoQuantitiesMatch = .false.
    if ( .not. DoQtysDescribeSameThing  ( a, b, options ) ) return
    ! Must we screen further?
    if ( .not. myStrict &
      & .and. &
      &.not. any( a%template%quantityType == screenedtypes ) &
      & ) then
      DoQuantitiesMatch = .true.
      return
    endif
    if ( .not. DoHGridsMatch ( a, b ) ) return
    if ( .not. DoVGridsMatch ( a, b ) ) return
    if ( .not. DoFGridsMatch ( a, b ) ) return
    DoQuantitiesMatch = .true.
  end function DoQuantitiesMatch

  ! -------------------------------------- DoVectorsMatch --------------
  logical function DoVectorsMatch ( a, b, verbose )
    ! Returns true if vectors are essentially the same in nature, even
    ! if names of quantities differ.
    type(Vector_T), intent(in) :: A ! First vector
    type(Vector_T), intent(in) :: B ! Second vector
    logical, optional, intent(in) :: verbose ! Say why not if they don't

    ! Local variables
    logical :: myVerbose
    integer :: Q                        ! Loop counter

    ! Exectuable code
    myVerbose = .false.
    if ( present(verbose) ) myVerbose = verbose
    DoVectorsMatch = .false.
    if ( a%template%noQuantities /= b%template%noQuantities ) then
      if ( myVerbose ) call outputnamedValue( 'noQuantities', &
        & (/ a%template%noQuantities, b%template%noQuantities /) )
      return
    endif
    do q = 1, a%template%noQuantities
      if ( .not. DoQuantitiesMatch ( &
        & a%quantities(q), b%quantities(q) ) ) then
        if ( myVerbose ) then
          call outputnamedValue( 'q', q )
          call outputnamedValue( 'doQuantitiesDesc', &
            & DoQtysDescribeSameThing  ( a%quantities(q), b%quantities(q), &
            & options='-d' ) )
          call outputnamedValue( 'doHGridsMatch', &
            & DoHGridsMatch ( a%quantities(q), b%quantities(q) ) )
          call outputnamedValue( 'doVGridsMatch', &
            & DoVGridsMatch ( a%quantities(q), b%quantities(q) ) )
          call outputnamedValue( 'doFGridsMatch', &
            & DoFGridsMatch ( a%quantities(q), b%quantities(q) ) )
          call dump( a%quantities(q) )
          call dump( b%quantities(q) )
        return
        endif
      endif
    end do
    DoVectorsMatch = .true.
  end function DoVectorsMatch

  ! ------------------------------------------  DoVGridsMatch_Vec  -----
  logical function DoVGridsMatch_Vec ( A, B, RelativeError, Precision )
    ! Returns true if quantities have same vGrid information
    use HyperSlabs, only: EssentiallyEqual
    use Dump_0, only: Dump
    type (vectorValue_T), intent(in) :: A ! First quantity
    type (vectorValue_T), intent(in) :: B ! Second quantity
    real(rv), optional, intent(in)   :: RelativeError ! May differ by this rel amount
    real(rv), optional, intent(in)   :: Precision ! May differ by this abs amount
    real, parameter                  :: defaultRelativeError = 1.e-9
    logical                          :: TestForSurfs
    ! logical, parameter               :: DEEBUG = .true.
    ! Executable code
    ! Can we do this?
    ! doVGridsMatch_Vec = a%template%VGridIndex /= 0 .and. &
    !                   & a%template%VGridIndex == b%template%VGridIndex
    ! if ( doVGridsMatch_Vec ) return
    doVGridsMatch_Vec = .false.
    if ( DEEBUG ) then
      call outputNamedValue ( 'noSurfs', (/ a%template%noSurfs, b%template%noSurfs /) )
      call outputNamedValue ( 'verticalCoordinate', (/ a%template%verticalCoordinate, b%template%verticalCoordinate /) )
      call outputNamedValue ( 'coherent', (/ a%template%coherent, b%template%coherent /) )
      call outputNamedValue ( 'regular', (/ a%template%regular, b%template%regular /) )
      call outputNamedValue ( 'noInstances', (/ a%template%noInstances, b%template%noInstances /) )
      call dump ( a%template%surfs, 'a%template%surfs' )
      call dump ( b%template%surfs, 'b%template%surfs' )
    endif
    if ( a%template%noSurfs /= b%template%noSurfs ) return
    if ( a%template%verticalCoordinate /= &
      &  b%template%verticalCoordinate ) return
    if ( a%template%coherent .neqv. b%template%coherent ) return
    if ( a%template%regular .neqv. b%template%regular ) return
    if ( ( .not. a%template%coherent) .and. &
      &  ( a%template%noInstances /= b%template%noInstances ) ) return
    ! Do we allow the corresponding surfaces any leeway? Relative or absolute?
    if ( present(RelativeError) ) then
      TestForSurfs = all ( abs(a%template%surfs - b%template%surfs) <= &
        & RelativeError*max(&
        & maxval( abs(a%template%surfs) ), maxval( abs(b%template%surfs ) ) &
        & ) &
        & )
    elseif ( present(Precision) ) then
      TestForSurfs = all ( abs(a%template%surfs - b%template%surfs) <= Precision )
    elseif ( defaultRelativeError > 0. ) then
      TestForSurfs = all ( abs(a%template%surfs - b%template%surfs) <= &
        & defaultRelativeError*max(&
        & maxval( abs(a%template%surfs) ), maxval( abs(b%template%surfs ) ) &
        & ) &
        & )
    else
      TestForSurfs = all (  essentiallyEqual ( a%template%surfs, &
                                      b%template%surfs ) )
    endif
    if ( .not. TestForSurfs ) return
    if (.not. a%template%regular ) then
      if ( any(a%template%surfIndex /= b%template%surfIndex) .or. &
        &  any(a%template%chanIndex /= b%template%chanIndex) ) return
    end if

    doVGridsMatch_Vec = .true.
  end function DoVGridsMatch_Vec

  ! --------------------------------------- DoFGridsMatch --------------
  logical function DoFGridsMatch ( a, b, sizeOnly )
    ! Returns true if the quantities have the same fGrid information
    type ( VectorValue_T ), intent(in) :: A ! First quantity
    type ( VectorValue_T ), intent(in) :: B ! Second quantity
    logical, intent(in), optional :: SIZEONLY

    ! Local parameters
    real (r8), parameter :: FTOL = 1.0e-3 ! 1 kHz
    logical :: MYSIZEONLY

    ! Executable code
    DoFGridsMatch = .false.
    mySizeOnly = .false.
    if ( present ( sizeOnly ) ) mySizeOnly = sizeOnly
    if ( a%template%noChans /= b%template%noChans ) return
    if ( .not. mySizeOnly ) then
      if ( a%template%frequencyCoordinate /= b%template%frequencyCoordinate ) return
      select case ( a%template%frequencyCoordinate )
      case ( l_none )
      case ( l_channel )
        if ( a%template%signal /= b%template%signal ) return
        if ( a%template%sideband /= b%template%sideband ) return
      case default
        if ( associated ( a%template%frequencies ) .neqv. &
          &  associated ( b%template%frequencies ) ) return
        if ( associated ( a%template%frequencies ) ) then
          if ( any ( shape(a%template%frequencies) /= &
            & shape(b%template%frequencies) ) ) return
          if ( any ( abs ( a%template%frequencies - &
            & b%template%frequencies ) > fTol ) ) return
        end if
      end select
    end if
    DoFGridsMatch = .true.
  end function DoFGridsMatch

  ! ---------------------------------- DoQtysDescribeSameThing ----
  logical function DoQtysDescribeSameThing ( a, b, options )
    ! Returns true if the quantities describe the same geophysical
    ! parameter, albeit at a different resolution perhaps
    
    ! Note:
    ! We have made this a lot more lenient than it was originally
    ! Unless you set options to include 's'[trict], it only checks 
    ! [type, logbasis, verticalcoord, signal, molecule)
    ! and then only if quantitype is one of eight types:
    ! [vmr, columnbundance, rhi, isotoperatio, radiance, 
    ! calsidebandfraction, limbsidebandfraction, or Tscat]
    ! Any other quantityTypes will be assumed to match if they
    ! are the same type, no further checks being required
    type ( VectorValue_T ), intent(in) :: A ! First quantity
    type ( VectorValue_T ), intent(in) :: B ! Second quantity
    character(len=*), optional, intent(in)      :: options ! e.g. 's'[trict]
    ! logical, optional, intent(in)      :: Strict ! Must every attribute match?
    !                                   Or just the ones we actually write out
    ! Local variables
    logical :: myStrict, debug
    ! logical, parameter :: DEEBUG = .true.
    integer, parameter :: NSCREENEDTYPES = 8
    integer, dimension(NSCREENEDTYPES) :: SCREENEDTYPES
    ! Executable
    screenedtypes = (/ &
      & l_vmr, l_columnabundance, l_rhi, l_isotoperatio, l_radiance, l_calsidebandfraction, l_limbsidebandfraction, l_Tscat &
      & /)
    myStrict = .false.
    if ( present(options) ) myStrict = index(options, 's') > 0
    debug = DEEBUG
    if ( present(options) ) debug = index(options, 'd') > 0

    DoQtysDescribeSameThing = .false.
    if ( debug ) then
      call outputNamedValue( 'myStrict', myStrict )
      call outputNamedValue( 'quantityType', (/ a%template%quantityType, b%template%quantityType /) )
      call outputNamedValue( 'logBasis', (/ a%template%logBasis, b%template%logBasis /) )
      call outputNamedValue( 'verticalCoordinate', (/ a%template%verticalCoordinate, b%template%verticalCoordinate /) )
      call outputNamedValue( 'unit', (/ a%template%unit, b%template%unit /) )
      call outputNamedValue( 'signal', (/ a%template%signal, b%template%signal /) )
      call outputNamedValue( 'molecule', (/ a%template%molecule, b%template%molecule /) )
      call outputNamedValue( 'sideband', (/ a%template%sideband, b%template%sideband /) )
      call outputNamedValue( 'instrumentModule', (/ a%template%instrumentModule, b%template%instrumentModule /) )
      call outputNamedValue( 'radiometer', (/ a%template%radiometer, b%template%radiometer /) )
    endif
    if ( a%template%quantityType /= b%template%quantityType ) return
    ! Must we screen further?
    if ( .not. myStrict &
      & .and. &
      &.not. any( a%template%quantityType == screenedtypes ) &
      & ) then
      DoQtysDescribeSameThing = .true.
      return
    endif
    if ( a%template%logBasis .neqv. b%template%logBasis ) return
    if ( a%template%verticalCoordinate /= b%template%verticalCoordinate ) return
    if ( a%template%signal /= b%template%signal ) return
    if ( a%template%molecule /= b%template%molecule ) return
    if ( myStrict ) then
      ! Because we chose not to write these as attributes
      ! E.g., in WriteVectorAsHDF5L2PC, or we neglected to
      ! read them in ReadOneVectorFromHDF5 we dare not automatically check
      if ( a%template%unit /= b%template%unit ) return
      if ( a%template%sideband /= b%template%sideband ) return
      if ( a%template%instrumentModule /= b%template%instrumentModule ) return
      if ( a%template%radiometer /= b%template%radiometer ) return
    endif
    DoQtysDescribeSameThing = .true.

  end function DoQtysDescribeSameThing

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ManipulateVectorQuantities
  
! $Log$
! Revision 2.54  2017/11/03 19:58:24  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.53  2016/06/03 02:33:06  vsnyder
! If both templates have pointers to their hGrids, verify that the hGrids
! have the same type.
!
! Revision 2.52  2015/10/28 00:36:39  vsnyder
! Confine WindowStart:WindowFinish with 1:noInstances
!
! Revision 2.51  2015/09/22 23:12:19  vsnyder
! Remove unused PHYQ_Angle USE reference
!
! Revision 2.50  2015/08/25 17:16:07  vsnyder
! In FindOneClosestInstance, determine whether to use the value or
! geolocation depending upon whether the type is PhiTan; eliminate the
! UseValue dummy argument.  In FindInstanceWindow, allow PhiWindow to be
! a tuple, with the first element giving the number of profiles/MAFs before
! the tangent point, and the second giving the number after.
!
! Revision 2.49  2015/07/29 00:27:28  vsnyder
! Convert Phi from pointer to allocated
!
! Revision 2.48  2015/06/19 00:11:01  pwagner
! Intercept and print if about to use unassociated phi in FindOneClosestInstance
!
! Revision 2.47  2015/02/05 21:43:08  vsnyder
! Don't use Phi for unstacked quantities
!
! Revision 2.46  2014/10/30 01:42:46  vsnyder
! Publish DoQuantitiesMatch
!
! Revision 2.45  2014/04/24 23:51:59  pwagner
! Depending on horizontalCoordinate of reference quantity, may FindOneClosestInstance in time, not phi
!
! Revision 2.44  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.43  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.42  2013/03/01 01:06:06  pwagner
! Get R8 from MLSKinds
!
! Revision 2.41  2012/07/19 03:33:18  vsnyder
! Pass nChan=quantity%template%noChans to CreateBlock
!
! Revision 2.40  2012/02/13 23:24:35  pwagner
! DoQuantitiesMatch takes options string; DoQtysDescribeSameThing more lenient
!
! Revision 2.39  2012/02/10 23:52:35  vsnyder
! Cannonball polishing
!
! Revision 2.38  2012/02/02 01:09:11  pwagner
! DoQuantitiesMatch now works properly when tested
!
! Revision 2.37  2011/08/29 21:29:41  pwagner
! Granted some leeway in matching vgrids
!
! Revision 2.36  2011/06/16 20:16:54  vsnyder
! Cannonball polishing
!
! Revision 2.35  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.34  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.33  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.32  2007/08/13 17:37:42  pwagner
! Push some procedures onto new MLSCallStack
!
! Revision 2.31  2006/03/03 23:05:50  pwagner
! Changed interface to AnyGoodDataInQty
!
! Revision 2.30  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.29  2004/12/13 20:28:44  vsnyder
! Changed DoVGridsMatch to a generic, with a DoVGridsMatch_Vec specific.
!
! Revision 2.28  2004/09/25 00:15:48  livesey
! Bug fixes etc. in FillWithCombinedChannels
!
! Revision 2.27  2004/09/24 17:55:41  livesey
! Gained ManipulateVectorQuantities from fill
!
! Revision 2.26  2004/01/24 01:01:48  livesey
! Improvements to DoFGridsMatch
!
! Revision 2.25  2004/01/23 05:37:01  livesey
! Added DoVectors/QuantitiesMatch
!
! Revision 2.24  2003/08/28 00:44:43  livesey
! Added sizeOnly option to DoFGridsMatch
!
! Revision 2.23  2003/07/07 20:21:34  livesey
! Now uses the FindClosestMatch function
!
! Revision 2.22  2003/01/26 04:42:20  livesey
! Added handling of profiles/angle units for phiWindow
!
! Revision 2.21  2002/11/22 01:07:13  vsnyder
! Delete USE'd but unreferenced symbols
!
! Revision 2.20  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.19  2002/09/19 00:30:36  pwagner
! Added AnyGoodDataInQty; set DoQtysDescribeSameThing at its start
!
! Revision 2.18  2002/09/10 20:47:44  livesey
! Added DoQtysDescribeSameThing
!
! Revision 2.17  2002/08/23 01:24:18  livesey
! Added DoFGridsMatch
!
! Revision 2.16  2002/07/25 08:43:19  mjf
! Initialised DoHGridsMatch to .false. at start of function.
!
! Revision 2.15  2002/07/17 06:01:27  livesey
! Fixed bugs in DoH/VGrids match
!
! Revision 2.14  2002/06/12 16:53:32  livesey
! Tidied up some public/private stuff
!
! Revision 2.13  2002/06/12 16:50:39  livesey
! Added findInstanceWindow
!
! Revision 2.12  2002/02/06 01:32:58  livesey
! Rewrote FindOneClosestInstance and FindClosestInstances to reflect the
! way in which they are mostly called, and to fix a bug.
!
! Revision 2.11  2001/11/08 01:05:06  livesey
! Fixed a minor sort of bug in FindOneClosestQuantity
!
! Revision 2.10  2001/09/14 18:02:52  livesey
! Bug fix in FindOneClosestInstance and FindClosestInstances.
! Will probably come back to these and rewrite them some time.
!
! Revision 2.9  2001/09/11 01:27:27  livesey
! Bug fixes
!
! Revision 2.8  2001/09/09 21:17:30  livesey
! Imported FindOneClosestInstance from branch
!
! Revision 2.7.2.1  2001/09/08 23:46:40  livesey
! Added FindOneClosestInstance
!
! Revision 2.7  2001/05/11 00:03:41  livesey
! Fixed but with DoVGridsMatch
!
! Revision 2.6  2001/05/10 23:29:27  livesey
! Added DoHGridsMatch and DoVGridsMatch
!
! Revision 2.5  2001/03/08 02:21:08  livesey
! Fixed bug, wasn't setting minloc!
!
! Revision 2.4  2001/03/02 01:31:36  livesey
! Regular commit
!
! Revision 2.3  2001/02/27 17:18:20  livesey
! Moved ValidateVectorQuantity into vectors module
!
