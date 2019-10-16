! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SubsetModule

  ! This module deals with all the user requests to 'subset' quantities
  ! through the mask field.  In the first instance, it was abstracted
  ! out of the retrieval module.

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  
  public :: RESTRICTRANGE, SETUPSUBSET, SETUPFLAGCLOUD
  public :: APPLYMASKTOQUANTITY, UPDATEMASK

  ! Local parameters
  character(len=32), private, parameter :: WRONGUNITS = &
    &                           'The wrong units were supplied'
    
contains ! ========= Public Procedures ============================

  ! ---------------------------------------------- ApplyMaskToQuantity -----
  ! Turn height node, or surfs node, and other tree nodes and flags
  ! into range of surfaces, instances, etc.
  ! Then possibly create and set qty mask over the *complement* of the range
  ! I.e., the range will be left unmasked so its values can be
  ! used, Filled, or whatever.

  ! We intend to make this the standard routine Subset and Fill commands
  ! call to restrict heights, channels, and horizontal instances
  ! instead of having each command process the l2cf fields
  ! separately
  subroutine ApplyMaskToQuantity ( qty, rad, ptan, opticaldepth, a, &
    & opticaldepthcutoff, maxvalue, minvalue, missingvalue, &
    & heightrange, whererange, &
    & ignore, reverse, additional, expandMask, reset, &
    & maskbit, heightnode, surfnode, instancesnode, channelsnode, got )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Dump_0, only: Dump
    use Expr_M, only: Expr, Getindexflagsfromlist
    use FillUtils_1, only: ByManipulation, M_All
    use Highoutput, only: OutputNamedValue
    use Init_Tables_Module, only: F_Height, &
      & F_MinValue, F_MaxValue, F_MissingValue, &
      & F_Ptanquantity, F_Quantity, F_Surface
    use Init_Tables_Module, only: L_None, L_Pressure, &
      & L_Zeta
    use Intrinsic, only: Phyq_Dimensionless, Phyq_Length, &
      & Phyq_Mifs, Phyq_Pressure
    use MLSKinds, only: R8, Rv
    use MLSStringLists, only: SwitchDetail
    use MLSFinds, only: FindFirst, FindLast
    use MLSStrings, only: TrueList
    use Output_M, only: Output
    use Toggles, only: Switches
    use Tree, only: Nsons, Subtree, Node_Id
    use Tree_Types, only: N_Colon_Less, N_Less_Colon, &
      & N_Less_Colon_Less
    use VectorsModule, only: VectorValue_T, &
      & ClearMask, CloneVectorQuantity, CreateMask, Dump, &
      & ReverseMask, SetMask
    ! Args
    type (VectorValue_T), pointer :: QTY
    type (VectorValue_T), pointer :: RAD
    type (VectorValue_T), pointer :: PTAN
    type (VectorValue_T), pointer :: OPTICALDEPTH
    type (VectorValue_T), pointer :: A
    real(r8), intent(in)          :: OPTICALDEPTHCUTOFF
    real(r8), intent(in)          :: MAXVALUE
    real(r8), intent(in)          :: MINVALUE
    real(r8), intent(in)          :: MISSINGVALUE
    character(len=*), intent(in)  :: HEIGHTRANGE
    integer, intent(in)           :: WHERERANGE
    logical, intent(in)           :: IGNORE
    logical, intent(in)           :: REVERSE
    logical, intent(in)           :: ADDITIONAL
    logical, intent(in)           :: ExpandMask
    logical, intent(in)           :: RESET
    integer, intent(in)           :: MASKBIT             ! Bits corresponding to Mask
    integer, intent(in)           :: HEIGHTNODE          ! Tree node
    integer, intent(in)           :: SURFNODE            ! Tree node
    integer, intent(in)           :: INSTANCESNODE       ! Tree node
    integer, intent(in)           :: CHANNELSNODE        ! Tree node
    logical, dimension(:), intent(in) :: Got

    ! Local variables
    type (VectorValue_T), pointer :: B
    logical, dimension(:), pointer :: CHANNELS ! Are we dealing with these channels
    integer :: CHANNEL
    integer :: COORDINATE
    logical :: DEEBUG
    logical :: DOTHISCHANNEL          ! Flag
    logical :: DOOKHEIGHTS            ! Flag
    logical :: DOTHISHEIGHT           ! Flag
    logical, dimension(:), pointer :: doThisInstance => null()
    integer :: FIRSTOK
    integer :: HEIGHT
    integer :: HEIGHTUNIT             ! Unit for heights command
    integer :: I
    integer :: IND
    integer :: INSTANCE
    integer :: INSTANCEOR1            ! For coherent quantities
    integer :: j
    integer :: LASTOK
    integer :: NODE
    character(len=128) :: NUMBERSLIST ! for a dump
    integer :: ODCUTOFFHEIGHT
    character(len=1), dimension(:), pointer :: ORIGINALMASK
    integer :: RANGEID
    integer, dimension(1) :: S1
    integer, dimension(1) :: S2
    integer               :: SS1
    integer               :: SS2
    integer :: SCANDIRECTION
    integer :: SON
    logical :: SPREADFLAG
    integer :: STATUS
    real(r8), dimension(:), pointer :: THESEHEIGHTS ! Subset of heights
    integer :: UNITS(2)               ! Units returned by expr
    real(r8) :: VALUE(2)              ! Value returned by expr
    logical :: VERBOSE
    type (VectorValue_T) :: WHEREAISOK ! where(a) /= 0
    ! Executable
    nullify ( channels, doThisInstance )
    verbose = ( switchDetail(switches,'subset') > -1 )
    DEEBUG = ( switchDetail(switches,'subset') > 0 )
    ! Do we have anything to do here?
    if ( all( (/ heightNode, surfNode, InstancesNode, channelsNode /) < 1) .and. &
      & .not. associated(opticalDepth) .and. .not. associated(a) .and. &
      & .not. any( (/reset, ignore/) ) ) then
      if ( verbose ) &
        & call output( 'No change needed to quantity mask', advance='yes' )
      return
    endif
    ! Process the instances field.
    call Allocate_test ( doThisInstance, qty%template%noInstances, &
      & 'doThisInstance', ModuleName )
    if ( DEEBUG ) then
      call outputNamedValue( 'MaskBit', MaskBit )
      call outputNamedValue( 'InstancesNode', InstancesNode )
      call outputNamedValue( 'HeightNode', HeightNode )
      call outputNamedValue( 'ChannelsNode', ChannelsNode )
      call outputNamedValue( 'SurfNode', SurfNode )
      call outputNamedValue( 'reset', reset )
      call outputNamedValue( 'ignore', ignore )
    end if
    doThisInstance = ( InstancesNode == 0 ) ! Unless specified, treat all instances alike
    if ( InstancesNode > 0 ) then
      call GetIndexFlagsFromList ( InstancesNode, doThisInstance, status )
      if ( status /= 0 ) call announceError ( 0, &
        & 'There was a problem with the instances field' )
      if ( DEEBUG ) then
        call trueList ( doThisInstance, numbersList )
        call outputNamedValue( 'Instances', trim(numbersList) )
        call dump( doThisInstance, name='Instances' )
      end if
    end if

    ! Process the channels field.
    if ( qty%template%frequencyCoordinate /= l_none .or. &
         associated(rad) ) then
      !??? Someday think about the low bound for channels using ???
      !??? lbound(spectrometerTypes(signals(...%template%signal)%spectrometerType)%Frequencies,1) ???
      if ( associated(rad) ) then
        call Allocate_test ( channels, rad%template%noChans, 'channels', &
        & ModuleName )
      else
        call Allocate_test ( channels, qty%template%noChans, 'channels', &
          & ModuleName )
      end if
      if ( channelsNode > 0 ) then     ! This subset is only for some channels
        call GetIndexFlagsFromList ( channelsNode, channels, status, &
          !??? lbound(channels,1) is always 1
          & lower=lbound(channels,1) )
        if ( status /= 0 ) call announceError ( channelsNode, &
          & 'There was a problem with the channels field' )
        if ( DEEBUG ) then
          call trueList ( channels, numbersList )
          call outputNamedValue( 'Channels', trim(numbersList) )
          call dump( channels, name='Channels' )
        end if
      else
        channels = .true.             ! Apply this to all channels
      end if
    end if

    ! Preprocess the height stuff.  
    heightUnit = phyq_dimensionless
    if ( heightNode > 0 ) then
      do j = 2, nsons(heightNode)
        call expr ( subtree(j,heightNode), units, value )
        ! Make sure the range has non-dimensionless units -- the type
        ! checker only verifies that they're consistent.  We need to
        ! check each range separately, because the units determine the
        ! scaling of the values.
        if ( all(units == phyq_dimensionless) ) call announceError ( &
          & subtree(j,heightNode), 'No height units', f_height )
        ! Check consistency of units -- all the same, or dimensionless. The
        ! type checker verifies the consistency of units of ranges, but not
        ! of array elements.
        do i = 1, 2
          if ( heightUnit == phyq_dimensionless ) then
            heightUnit = units(i)
          else if ( units(i) /= phyq_dimensionless .and. &
            &       units(i) /= heightUnit ) then
            call announceError ( heightNode, 'Inconsistent height units', f_height )
          end if
        end do
      end do
      ! Check for correct units
      if ( heightUnit == phyq_pressure ) then
        if ( qty%template%minorFrame .and. .not. associated(ptan) ) &
          & call announceError ( heightNode, 'Needed ptan', f_height, f_ptanQuantity )
      else if ( heightUnit /= phyq_length ) then
        call announceError ( heightNode, 'minor frame height units must be MIFs', f_height )
      end if
    end if

    ! Create the mask if it doesn't exist
    if ( .not. associated( qty%mask ) ) call CreateMask ( qty )

    ! Make a space to save the original values if doing an additional mask
    if ( additional ) then
      nullify ( originalMask )
      call Allocate_test ( originalMask, qty%template%instanceLen, &
        & 'originalMask', ModuleName )
    end if

    if ( associated(a) ) then
      nullify( b )
      call CloneVectorQuantity( whereAIsOK, qty )
      spreadFlag = qty%template%NoInstances /= a%template%NoInstances
      ! Evaluate where(a)
      call ByManipulation ( whereAIsOK, a, b, &
        & whereRange, key=0, ignoreTemplate=.true., &
        & spreadflag=spreadFlag, dimList=' ', &
        & c=0._rv )
      if ( verbose ) then
        call dump( whereAIsOK, details=1 )
      endif
    endif
    ! Now we loop over the instances
    do instance = 1, qty%template%noInstances
      ! Have we been asked to expand the original mask?
      if ( expandMask ) then
        where (ichar(qty%mask ( :, instance )) /= 0 )
          qty%mask ( :, instance ) = char(M_All)
        end where
        cycle
      endif
      ! Possibly save original mask
      if ( additional ) originalMask = qty%mask ( :, instance )
      if ( .not. doThisInstance(instance) ) then
        ! Because this instance is outside the instances=.. field 
        ! we will mask this instance (because subset picks out
        ! a portion of a quantity not to mask)
        !??? Make sure mask bit numbers begin at 1, even when
        !??? channel numbers don't.
        call SetMask ( qty%mask(:,instance), &
          & what=maskBit )
        ! If this is supposed to be an 'additional' mask, merge in the
        ! original value
        if ( additional ) &
          & qty%mask(:,instance) = char ( ior ( &
          & ichar ( qty%mask(:,instance) ), ichar ( originalMask ) ) )
        cycle ! instance
      end if

      if ( associated(a) ) then
        call SetMask ( qty%mask(:,instance), &
          & what=maskBit )
        ! Find first, last heights at which whereAIsOK is 1
        firstOK = FindFirst( whereAIsOK%values(:, instance), 1._rv )
        lastOK  = FindLast ( whereAIsOK%values(:, instance), 1._rv )
        if ( verbose ) then
          call outputnamedValue( 'firstOK', firstOK )
          call outputnamedValue( 'firstOK', firstOK )
          call outputnamedValue( 'heightRange', heightRange )
        endif
        ! Now our interpretation:
        ! If heightRange is missing, we will unmask only heights
        ! which are OK
        ! If 'above', we will unmask all heights above the lowest OK
        ! If 'below', we will unmask all heights below the highest OK
        ! These height ranges are set by ss1 and ss2
        select case (heightRange)
        case ('above')
          ss1 = firstOK
          ss2 = qty%template%noSurfs
          doOkHeights = .false.
        case ('below')
          ss1 = 1
          ss2 = lastOK
          doOkHeights = .false.
        case default
          ss1 = 1
          ss2 = qty%template%noSurfs
          doOkHeights = .true.
        end select
        if ( verbose ) then
          call outputnamedValue( 'ss1', ss1 )
          call outputnamedValue( 'ss2', ss2 )
        endif
        ! If ss1 or ss2 is 0, bail; otherwise unset mask bits
        if ( ss1 > 0 .and. ss2 > 0 ) then
          do channel = 1, qty%template%noChans
            doThisChannel = .true.
            if ( associated(channels) ) doThisChannel = channels(channel)
            if ( doThisChannel ) then
              do height = ss1, ss2
                !??? Make sure mask bit numbers begin at 1, even when
                !??? channel numbers don't.
                if ( doOkHeights ) then
                  if ( whereAIsOK%values( &
                    & channel+qty%template%noChans*(height-1), instance&
                    & ) /= 1._rv ) cycle
                endif
                call ClearMask ( qty%mask(:,instance), &
                  & (/ channel+qty%template%noChans*(height-1) /), &
                  & what=maskBit )
              end do                  ! Height loop
            end if                    ! Do this channel
          end do                      ! Channel loop
        endif
        ! If this is supposed to be an 'additional' mask, merge in the
        ! original value
        if ( additional ) &
          & qty%mask(:,instance) = char ( ior ( &
          & ichar ( qty%mask(:,instance) ), ichar ( originalMask ) ) )
        cycle ! instance
      end if

      if ( associated(rad) ) then
        do height = 1, qty%template%noSurfs
          ! Unmask qty where any specified channel in rad is unmasked
          qty%mask(height,instance) = char(255)
          do channel = 1, rad%template%noChans
            if ( ignore .neqv. channels(channel) ) &
              & qty%mask(height,instance) = char( &
                & iand( ichar(qty%mask(height,instance)), &
                      & ichar(rad%mask(channel+rad%template%noChans*(height-1),instance)) ) )
          end do
        end do
        ! If this is supposed to be an 'additional' mask, merge in the
        ! original value
        if ( additional ) &
          & qty%mask(:,instance) = char ( ior ( &
          & ichar ( qty%mask(:,instance) ), ichar ( originalMask ) ) )
        cycle ! instance
      end if

      instanceOr1 = instance

      if ( qty%template%coherent ) then
        theseHeights => qty%template%surfs(:,1)
        coordinate = qty%template%verticalCoordinate
        instanceOr1 = 1
      else if ( qty%template%minorFrame .and. heightUnit == phyq_pressure ) then
        if ( ptan%template%instrumentModule /= qty%template%instrumentModule ) &
          & call AnnounceError ( 0, 'ptan is for wrong module', f_ptanQuantity, &
          & f_quantity )
        theseHeights => ptan%values(:,instance)
        coordinate = l_zeta
      else
        theseHeights => qty%template%surfs(:,instance)
        coordinate = qty%template%verticalCoordinate
      end if

      ! Now, make sure for the channels we're considering that the
      ! default is to ignore all, unless we've not got a heights or ignore
      ! in which case we default to use all
      ! if ( any ( (/ got(f_height), got(f_surface), ignore /) ) ) then
      if ( ignore .or. heightNode > 0 .or. surfNode > 0 ) then
        do channel = 1, qty%template%noChans
          doThisChannel = .true.
          if ( associated(channels) ) doThisChannel = channels(channel)
          if ( doThisChannel ) then
            do height = 1, qty%template%noSurfs
              !??? Make sure mask bit numbers begin at 1, even when
              !??? channel numbers don't.
              call SetMask ( qty%mask(:,instance), &
                & (/ channel+qty%template%noChans*(height-1) /), &
                & what=maskBit )
            end do                    ! Height loop
          end if                      ! Do this channel
        end do                        ! Channel loop
      else
        ! If not got heights and/or ignore, default must be 
        if ( reset ) then
          do channel = 1, qty%template%noChans
            doThisChannel = .true.
            if ( associated(channels) ) doThisChannel = channels(channel)
            if ( doThisChannel ) then
              do height = 1, qty%template%noSurfs
                !??? Make sure mask bit numbers begin at 1, even when
                !??? channel numbers don't.
                call ClearMask ( qty%mask(:,instance), &
                  & (/ channel+qty%template%noChans*(height-1) /), &
                  & what=maskBit )
              end do                  ! Height loop
            end if                    ! Do this channel
          end do                      ! Channel loop
        end if                        ! Reset specified
      end if                          ! Got heights or ignore

      ! Now go and `unmask' the ones we want to consider.  For coherent
      ! quantities we can simply loop over the indices of each height range.
      ! For incoherent ones we have to go through and mark each point
      ! appropriately.
      if ( heightNode > 0 .or. surfNode > 0 ) then
        if ( heightNode > 0 ) then
          node = heightNode
        else
          node = surfNode
        end if
        do j = 2, nsons(node)
          son = subtree ( j, node )
          rangeId = node_id ( son )
          if ( heightNode > 0 ) then
            ! Get values for this range
            call expr ( son, units, value )
            if ( DEEBUG ) call outputNamedValue( 'value', value(1) )
            ! Now maybe do something nasty to value to get in right units.
            if ( coordinate == l_zeta .and. heightUnit == phyq_pressure ) then
              value = -log10(value)
            else
              if ( coordinate /= qty%template%verticalCoordinate ) &
                & call AnnounceError ( son, 'vert. coord. wrong for template', f_height )
            end if
            if ( DEEBUG ) call outputNamedValue( '-log10(value)', value(1) )

            ! Do special things for coherent quantities
            if ( qty%template%coherent ) then
              if ( coordinate == l_pressure ) then
                s1 = minloc ( abs ( -log10(theseHeights) + log10(value(1)) ) )
                s2 = minloc ( abs ( -log10(theseHeights) + log10(value(2)) ) )
              else
                s1 = minloc ( abs ( theseHeights - value(1) ) )
                s2 = minloc ( abs ( theseHeights - value(2) ) )
              end if
              ! Did we try heightRange mechanism?
              select case ( heightRange )
              case ( 'above' )
                if ( theseHeights(s1(1)) < value(1) ) s1 = s1 + 1
              case ( 'below' )
                if ( theseHeights(s2(1)) > value(2) ) s2 = s2 - 1
              end select
            end if
          else
            ! Subset by surface (i.e. MIF) instead, much easier
            call expr ( son, units, value )
            if ( .not. all ( units == phyq_dimensionless .or. units == phyq_mifs ) ) &
              & call AnnounceError ( son, 'surface must be MIFs', f_surface )
            s1(1) = max ( min ( value(1), real(qty%template%noSurfs, r8) ), 1._r8 )
            s2(1) = max ( min ( value(2), real(qty%template%noSurfs, r8) ), 1._r8 )
          end if
          ! Now consider the open range issue
          select case ( rangeId )
          case ( n_colon_less )
            s1 = min ( s1 + 1, qty%template%noSurfs )
          case ( n_less_colon )
            s2 = max ( s2 - 1, 1 )
          case ( n_less_colon_less )
            s1 = min ( s1 + 1, qty%template%noSurfs )
            s2 = max ( s2 - 1, 1 )
          end select
          ! Did we try heightRange mechanism?
          select case ( heightRange )
          case ( 'above' )
            s2 = qty%template%noSurfs
          case ( 'below' )
            s2 = s1
            s1 = 1
          end select
          if ( DEEBUG ) then
            call outputNamedValue( 's1', s1(1) )
            call outputNamedValue( 's2', s2(1) )
          endif
          
          scanDirection = 0
          do channel = 1, qty%template%noChans
            doThisChannel = .true.
            if ( associated(channels) ) doThisChannel = channels(channel)
            if ( doThisChannel ) then

              ! Think about optical depth.  We want the cutoff to apply once
              ! and definitively, not have channels come in and out of use as
              ! a function of tangent height.  To do this we need to work out
              ! where we 'first' go optically thick. To do this we need to
              ! work out whether we're scanning up or down.
              if ( associated ( opticalDepth ) ) then
                ! Don't bother deducing direction if we alredy know
                if ( scanDirection == 0 ) then
                  ind = qty%template%noChans + channel
                  do height = 2, qty%template%noSurfs
                    scanDirection = scanDirection + merge ( 1, -1, &
                      & opticalDepth%template%surfs(height,instanceOr1) > &
                      & opticalDepth%template%surfs(height-1,instanceOr1) )
                  end do
                  ! Now convert it to +/-1.
                  if ( scanDirection == 0 ) scanDirection = 1 ! Default upscan
                  scanDirection = scanDirection / abs(scanDirection) 
                end if
                ! Now find `first' optically thick radiance look down from top
                if ( scanDirection == 1 ) then ! Scanning up
                  ind = qty%template%noChans*(qty%template%noSurfs-1) + channel
                  do odCutoffHeight = qty%template%noSurfs, 1, - 1 
                    if ( opticalDepth%values ( ind, instance ) > &
                      & opticalDepthCutoff ) exit
                    ind = ind - qty%template%noChans
                  end do
                else ! else scanning down
                  ind = channel ! qty%template%noChans + channel
                  do odCutoffHeight = 1, qty%template%noSurfs
                    if ( opticalDepth%values ( ind, instance ) > &
                      & opticalDepthCutoff ) exit
                    ind = ind + qty%template%noChans
                  end do
                end if
              end if

              if ( qty%template%coherent .or. surfNode > 0 ) then
                ! For coherent quantities and/or specified surfaces simply do a loop over a range.
                do height = s1(1), s2(1)
                  !??? Make sure mask bit numbers begin at 1, even when
                  !??? channel numbers don't.
                  ind = channel + qty%template%noChans*(height-1)
                  doThisHeight = .true.
                  if ( got(f_minValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) > minValue
                  if ( got(f_maxValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) < maxValue
                  if ( got(f_missingValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) == missingValue
                  if ( doThisHeight ) then
                    if ( associated( opticalDepth ) ) then
                      if ( scanDirection * ( height - odCutoffHeight ) > 0 ) &
                        & call ClearMask ( qty%mask(:,instance), (/ind/), what=maskBit )
                    else
                      call ClearMask ( qty%mask(:,instance), (/ind/), what=maskBit )
                    end if
                  end if
                end do                ! Height loop
              else
                ! For Incoherent quantities check each height individually
                ! We might as well consider open and non open ranges, but
                ! it's really rather unnecessary, as the chance of a ptan
                ! being exactly equal to the range given is remote, and the
                ! difference probably doesn't matter anyway.
                do height = 1, qty%template%noSurfs
                  ind = channel + qty%template%noChans*(height-1)
                  doThisHeight = .true.
                  if (any(rangeID==(/ n_less_colon,n_less_colon_less /))) then
                    doThisHeight = doThisHeight .and. theseHeights(height) > value(1)
                  else
                    doThisHeight = doThisHeight .and. theseHeights(height) >= value(1)
                  end if
                  if (any(rangeID==(/ n_colon_less,n_less_colon_less /))) then
                    doThisHeight = doThisHeight .and. theseHeights(height) < value(2)
                  else
                    doThisHeight = doThisHeight .and. theseHeights(height) <= value(2)
                  end if
                  if ( associated ( opticalDepth ) ) then
                    doThisHeight = doThisHeight .and. &
                      & scanDirection * ( height - odCutoffHeight ) > 0
                  end if
                  if ( got(f_minValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) > minValue
                  if ( got(f_maxValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) < maxValue
                  if ( got(f_missingValue) ) doThisHeight = doThisHeight .and. &
                    & qty%values( ind, instance ) == missingValue
                  if ( doThisHeight ) call ClearMask ( qty%mask(:,instance), &
                    & (/ ind /), what=maskBit )
                end do                ! Height loop
              end if                  ! Coherent
            end if                    ! Do this channel
          end do                      ! Channel loop
        end do                        ! Height entries in l2cf
      end if                          ! Got a height entry

      ! If the reverse flag is set, reverse the mask
      if ( reverse ) &
        & call ReverseMask( qty%mask(:,instance), what=maskBit )

      ! If this is supposed to be an 'additional' mask, merge in the
      ! original value
      if ( additional ) &
        & qty%mask(:,instance) = char ( ior ( &
        & ichar ( qty%mask(:,instance) ), ichar ( originalMask ) ) )
    end do                            ! Instance loop

    ! Tidy up
    call Deallocate_test ( channels, 'channels', ModuleName )
    call Deallocate_test ( dothisInstance, 'dothisInstance', ModuleName )
    if ( additional ) call Deallocate_test ( originalMask, &
      & 'originalMask', ModuleName )
          
  end subroutine ApplyMaskToQuantity

  ! -------------------------------------------------- RestrictRange ---
  ! This looks like it might be a nice operation. 
  ! However it comes without documentation, so who knows
  ! what it does or whether it works properly?

  ! So be warned--
  ! d o c u m e n t   y o u r   c o d e
  subroutine RestrictRange ( key, vectors )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_M, only: Expr
    use VectorsModule, only: M_Linalg, Vector_T, Vectorvalue_T, &
      & CreateMask, GetVectorQtyByTemplateIndex, GetVectorQuantityByType, &
      & SetMask
    use Init_Tables_Module, only: F_Quantity, F_Ptanquantity, F_Basisfraction, &
      & F_Minchannels, F_Signals, F_Measurements, F_Mask
    use Init_Tables_Module, only: L_Zeta, L_Radiance
    use ManipulateVectorQuantities, only: FindOneClosestInstance
    use MLSKinds, only: R8
    use MLSNumerics, only: Hunt
    use MLSSignals_M, only: Signal_T
    use MoreTree, only: Get_Field_Id
    use Parse_Signal_M, only: Expand_Signal_List
    use Tree, only: Nsons, Subtree, Decoration

    ! Dummy arguments
    integer, intent(in) :: KEY          ! Tree node
    type (Vector_T), dimension(:), intent(inout), target :: VECTORS

    ! Local variables

    integer :: CHANNEL                  ! Loop counter / channel index
    integer :: FIELD                    ! ID for field in l2cf line
    integer :: GSON                     ! Tree node
    integer :: INDEX                    ! Loop counter
    integer :: INDEX0                   ! Array index
    integer :: INDEX1                   ! Array index
    integer :: INSTANCE                 ! Instance index
    integer :: J                        ! Loop counter
    integer :: MAF                      ! Major frame index
    integer :: MASK                     ! Bits to mask
    integer :: MEASQTY                  ! Loop counter
    integer :: MINCHANNELS              ! Value of argument of same name
    integer :: NOGOODMIFS               ! How many MIFs meet our criteria for one chan
    integer :: NOINSTANCES              ! Number of instances
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOMEASUREMENTS           ! Number of relevant measurements
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOVALIDCHANNELS          ! How many channels meet our criteria
    integer :: QUANTITYINDEX            ! Index in database, not vector
    integer :: SON                      ! Tree node
    integer :: SURF                     ! Loop counter
    integer :: VECTORINDEX              ! Vector index for quantity

    integer, dimension(2) :: EXPRUNIT   ! Units for expression
    integer, dimension(:), pointer :: MAFSFORINSTANCE ! Which maf relevant for each inst.
    integer, dimension(:), pointer :: MEASQTYINDS ! Indices
    integer, dimension(:), pointer :: MIFPOINTINGS ! Which basis for each mif?

    logical :: ERRORFLAG                ! Set if problem
    logical :: FOUNDONE                 ! Flag for instance identification

    character(len=1), dimension(:), pointer :: MIFMASKS ! Part of rad%mask

    real(r8) :: BASISFRACTION           ! Value of argument of same name
    real(r8), dimension(2) :: EXPRVALUE ! Vaue for expression
    real(r8), dimension(:), pointer :: SURFS ! Surfaces for quantity

    type (Signal_T), dimension(:), pointer :: SIGNALS
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity
    type (VectorValue_T), pointer :: QTY ! The quantity to subset
    type (VectorValue_T), pointer :: RAD ! A radiance quantity
    type (Vector_T), pointer :: MEASUREMENTS ! Measurement vector

    ! Executable code
    nullify ( signals )
    ! Setup defaults
    basisFraction = 0.5
    minChannels = 1
    mask = m_linAlg
    ! Loop over arguments to this l2cf command
    do j = 2, nsons ( key )
      son = subtree ( j, key )
      field = get_field_id ( son )
      if ( nsons ( son )  > 1 ) gson = subtree ( 2, son )
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )
      case ( f_ptanQuantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )
      case ( f_minChannels )
        call expr ( gson, exprUnit, exprValue )
        minChannels = nint ( exprValue(1) )
      case ( f_basisFraction )
        call expr ( gson, exprUnit, exprValue )
        basisFraction = exprValue(1)
      case ( f_measurements )
        vectorIndex = decoration(decoration(gson))
        measurements => vectors ( vectorIndex )
      case ( f_signals )
        call Expand_Signal_List ( son, signals, errorFlag )
        if ( errorFlag ) call AnnounceError ( key, 'Bad signal spec.' )
      case ( f_mask )
        mask = GetMaskBit ( son )
      end select
    end do
    ! The parser has ensured that we 'got' all the fields we required,
    ! and we supplied defaults for the rest, so we're ready to go now.
    ! Some quick error checking
    if ( qty%template%verticalCoordinate /= l_zeta ) &
      & call AnnounceError ( key, 'Quantity is not on log pressure surfaces' )
    if ( .not. associated ( qty%mask ) ) call CreateMask ( qty )

    ! Identify all the measurement quantities we'll be using
    noMeasurements = size(signals)
    nullify ( measQtyInds )
    call Allocate_test ( measQtyInds, noMeasurements, 'measQtyInds', ModuleName )
    do measQty = 1, noMeasurements
      rad => GetVectorQuantityByType ( measurements, quantityType=l_radiance, &
        & signal=signals(measQty)%index, sideband=signals(measQty)%sideband )
      if ( rad%template%instrumentModule /= ptan%template%instrumentModule ) &
        & call AnnounceError ( key, 'Ptan/radiance not from the same module' )
      measQtyInds(measQty) = rad%index
    end do

    ! Now identify the relevant MAFs for each profile of our quantity
    ! This kind of involves running FindOneClosestInstance backwards
    noInstances = qty%template%noInstances
    noMAFS = ptan%template%noInstances
    nullify ( mafsForInstance )
    call Allocate_test ( mafsForInstance, noInstances, 'mafsForInstance', ModuleName )
    mafsForInstance = 0
    do maf = 1, noMAFs
      instance = FindOneClosestInstance ( qty, ptan, maf )
      mafsForInstance ( instance ) = maf
    enddo
    ! Deal with any we missed (probably at the ends)
    foundOne = .false.
    do instance = 1, noInstances
      if ( mafsForInstance ( instance ) /= 0 ) then
        foundOne = .true.
      else
        if ( foundOne ) then
          ! OK, make this gap contain the previous value,
          ! we know instance>1 here.
          mafsForInstance ( instance ) = mafsForInstance ( instance - 1 )
        else
          ! OK, still not seen on, make this use the first maf
          mafsForInstance ( instance ) = 1
        end if
      end if
    end do

    ! Setup some more arrrays
    nullify ( mifPointings )
    call Allocate_test ( mifPointings, ptan%template%noSurfs, &
      & 'mifPointings', ModuleName )

    ! Now we get down to the main work, we'll go through this instance by instance
    do instance = 1, noInstances
      maf = MAFsForInstance ( instance )
      if ( qty%template%coherent ) then
        surfs => qty%template%surfs ( :, 1 )
      else
        surfs => qty%template%surfs ( :, instance )
      end if
      ! Work out where each MIF points
      call Hunt ( surfs, ptan%values(:,maf), mifPointings, allowTopValue=.true. )
      ! Go through each surface/channel in the quantity and work out whether
      ! We're going to use it.
      surfLoop: do surf = 1, qty%template%noSurfs
        index0 = ( surf-1 ) * qty%template%noChans + 1
        index1 = index0 + qty%template%noChans - 1
        ! If were retrieving data for any 'channel' at this surface
        ! think about changing our minds, otherwise don't bother.
        if ( any ( iand ( ichar(qty%mask(index0:index1,instance)), &
          & m_linAlg ) == 0 ) ) then
          ! How many MIFs total above this surface?
          noMIFs = count ( mifPointings == surf )
          ! Now, how many channels have enough good channels (i.e. enough channels
          ! with more valid MIFs in this basis range than goodMIFs*basisFraction)
          noValidChannels = 0
          do measQty = 1, noMeasurements
            rad => measurements%quantities ( measQtyInds ( measQty ) )
            if ( associated ( rad%mask ) ) then
              channelLoop: do channel = 1, rad%template%noChans
                if ( .not. signals(measQty)%channels(channel) ) cycle channelLoop
                ! Point to the relevant masks for this radiance channel only
                mifMasks => rad%mask ( channel : rad%template%instanceLen : &
                  & rad%template%noChans, instance )
                noGoodMIFs = count ( mifPointings == surf .and. &
                  & iand( ichar(mifMasks), m_linAlg ) == 0 )
                if ( (1.0*noGoodMIFs)/noMIFs > basisFraction ) &
                  & noValidChannels = noValidChannels + 1
              end do channelLoop
            else
              noValidChannels = noValidChannels + 1
            end if
          end do ! Signal loop
          if ( noValidChannels < minChannels ) then
            ! We don't have enough information to retrieve this value,
            ! So mark it as not to be retrieved (or whatever bits have been supplied)
            do index = index0, index1
              call SetMask ( qty%mask(:,instance), (/ index /), what=mask )
            end do
          else
            ! We can retrieve this surface.  Not only that, I want this code
            ! to only consider the 'bottom' limit case, above that I don't
            ! want it to bother, so I'm going to escape from the surface loop here
            exit surfLoop
          end if                        ! Is this surface ok?
        end if                          ! Was planning to retrieve here
      end do surfLoop                   ! Surface loop
    end do                              ! Instance loop

    call Deallocate_test ( mifPointings, 'mifPointings', ModuleName )
    call Deallocate_test ( MAFsForInstance, 'MAFsForInstance', ModuleName )
    call Deallocate_test ( measQtyInds, 'measQtyInds', ModuleName )
  end subroutine RestrictRange

  ! -------------------------------------------------- SetupSubset ---
  ! (Documented in wiki)
  ! Called in response to Subset command. Process all its fields and
  ! set quantity mask.
  ! Mask setting business moved from here to workhorse subroutine.
  subroutine SetupSubset ( key, vectors )

    use Expr_M, only: Expr
    use Init_Tables_Module, only: Field_First, Field_Last
    use Init_Tables_Module, only: F_A, F_Additional, F_Channels, F_ExpandMask, &
      & F_Height, F_HeightRange, &
      & F_Ignore, F_Instances, &
      & F_Mask, F_MaxValue, F_MinValue, F_MissingValue, F_OpticalDepth, &
      & F_OpticalDepthCutoff, F_PTanQuantity, F_Quantity, F_RadianceQuantity, &
      & F_Reset, F_ResetAll, F_Reverse, F_SourceQuantity, F_Surface, F_Where
    use Init_Tables_Module, only: L_OpticalDepth, &
      & L_Radiance
    use Intrinsic, only: Phyq_Dimensionless, Phyq_Invalid
    use MLSKinds, only: R8
    use MLSStringlists, only: SwitchDetail
    use MoreTree, only: Get_Field_Id, Get_Boolean
    use Output_M, only: Output
    use String_Table, only: Display_String, Get_String
    use Toggles, only: Switches
    use Tree, only: Nsons, Sub_Rosa, Subtree, Decoration
    use VectorsModule, only: M_Linalg, Vector_T, VectorValue_T, &
      & CreateMask, DestroyVectorQuantityMask, DumpMask, &
      & GetVectorQtyByTemplateIndex

    integer, intent(in) :: KEY        ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    integer :: CHANNELSNODE           ! Tree node for channels values
    integer :: FIELD                  ! Field type from tree
    integer :: GSON                   ! Tree node
    integer :: HEIGHTNODE             ! Tree node for height values
    character(len=8) :: HEIGHTRANGE   ! 'above', 'below', or ' '
    integer :: J                      ! Subscripts, loop inductors
    integer :: INSTANCESNODE          ! Tree node for instance values
    integer :: MAINVECTORINDEX        ! Vector index of quantity to subset
    integer :: MANIPULATION
    integer :: MASKBIT                ! Bits corresponding to Mask
    integer :: MAXUNIT                ! Units for maxValue
    integer :: MINUNIT                ! Units for minValue
    integer :: QUANTITYINDEX          ! Index
    integer :: SON                    ! Tree node
    integer :: SURFNODE               ! Tree node
    integer :: TESTUNIT               ! Either vector%globalUnit or qty tmplt unit
    integer :: UNITS(2)               ! Units returned by expr
    integer :: VECTORINDEX            ! Index
    integer :: WHERERANGE             ! E.g., 'a > 100'

    real(r8) :: VALUE(2)              ! Value returned by expr
    real(r8) :: OPTICALDEPTHCUTOFF    ! Maximum value of optical depth to allow
    real(r8) :: MAXVALUE              ! Cutoff ranges
    real(r8) :: MINVALUE              ! Cutoff ranges
    real(r8) :: MISSINGVALUE              ! Cutoff ranges
    type (VectorValue_T), pointer :: A  ! The a of the 'where' field if needed
    type (VectorValue_T), pointer :: OPTICALDEPTH ! The opticalDepth quantity if needed
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity if needed
    type (VectorValue_T), pointer :: QTY  ! The quantity to mask
    type (VectorValue_T), pointer :: RAD  ! The quantity to mask
    type (VectorValue_T), pointer :: SOURCE  ! The sourceQty of mask if needed
    logical :: Got(field_first:field_last)   ! "Got this field already"
    logical :: IGNORE                 ! Flag
    logical :: RESET                  ! Flag
    logical :: RESETAll               ! Flag
    logical :: ADDITIONAL             ! Flag
    logical :: ExpandMask             ! Flag
    logical :: REVERSE                ! Flag

    ! Executable code
    nullify ( a, qty, ptan, opticalDepth, rad, source )
    got = .false.
    instancesNode = 0
    heightNode = 0
    heightRange = ' '
    surfNode = 0
    channelsNode = 0
    maxValue = -999.99
    minValue = -999.99
    missingValue = -999.99
    opticalDepthCutoff = -999.99
    ignore = .false.
    reset = .false.
    resetAll = .false.
    additional = .false.
    ExpandMask = .false.
    reverse = .false.
    maskBit = m_linalg
    minUnit = 0
    maxUnit = 0
    whereRange = 0
    do j = 2, nsons(key) ! fields of the "subset" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_a )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        a => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_quantity )
        mainVectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex(vectors(mainVectorIndex), quantityIndex)
      case ( f_ptanquantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_radiancequantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        rad => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_sourcequantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        source => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_channels )
        channelsNode = son
      case ( f_height )
        heightNode = son
      case ( f_heightRange )
        manipulation = sub_rosa ( gson )
        heightRange = ' '
        ! If heightRange field was present, it should have been one of
        ! 'a[bove]' meaning fill heights above supplied value
        ! 'b[elow]' meaning fill heights below supplied value
        call get_string ( manipulation, heightRange, strip=.true. )
        select case ( heightRange(1:1) )
        case ( 'a' )
          heightRange = 'above'
        case ( 'b' )
          heightRange = 'below'
        case ( ' ' )
          heightRange = ' '
        case default
          call AnnounceError ( key, &
            & 'invalid heightRange: ' // trim(heightRange) )
        end select
      case ( f_instances )
        instancesNode = son
      case ( f_mask )
        if ( .not. got(f_mask) ) maskBit = 0 ! clear default first time
        maskBit = GetMaskBit ( son )
      case ( f_opticalDepth )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        opticalDepth => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_opticalDepthCutoff )
        call expr ( subtree (2, son), units, value )
        if ( units(1) /= phyq_dimensionless ) &
          & call announceError ( subtree (2, son), WrongUnits, &
          & f_opticalDepthCutoff )
        opticalDepthCutoff = value(1)
      case ( f_maxValue )
        call expr ( subtree (2, son), units, value )
        maxUnit = units(1)
        maxValue = value(1)
      case ( f_minValue )
        call expr ( subtree (2, son), units, value )
        minUnit = units(1)
        minValue = value(1)
      case ( f_missingValue )
        call expr ( subtree (2, son), units, value )
        missingValue = value(1)
      case ( f_ignore )
        ignore = Get_Boolean ( son )
      case ( f_reset )
        reset = Get_Boolean ( son )
      case ( f_resetAll )
        resetAll = Get_Boolean ( son )
      case ( f_surface )
        surfNode = son
      case ( f_additional )
        additional = Get_Boolean ( son )
      case ( f_ExpandMask )
        ExpandMask = Get_Boolean ( son )
      case ( f_reverse )
        reverse = Get_Boolean ( son )
      case ( f_where )
        whereRange = sub_rosa ( gson )
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! SourceQuantity:
    ! Just copy the mask, if any, from the source
    if ( got(f_sourceQuantity) ) then
      if ( .not. associated( qty%mask ) ) call CreateMask ( qty )
      if ( .not. associated( source%mask ) ) then
        call destroyVectorQuantityMask ( qty )
        return
      endif
      qty%mask = source%mask
      return
    endif
    if ( got(f_radianceQuantity) ) then
      if ( rad%template%quantityType /= l_radiance ) &
        & call AnnounceError ( key, &
        & 'RadianceQuantity field is not a radiance quantity' )
      if ( any(got((/ f_height, f_mask, f_maxvalue, &
        & f_minvalue, f_opticaldepth, f_opticaldepthcutoff, f_ptanquantity, &
        & f_reset, f_reverse, f_surface /)))) call AnnounceError ( key, &
          & 'Subset from radiance allows only ' // &
          & 'additional, channel, ignore, and instances fields' )
      if ( .not. qty%template%minorFrame ) call AnnounceError ( key, &
        & 'Can only mask minor frame quantities from radiance' )
      if ( qty%template%noChans > 1 ) call AnnounceError ( key, &
        & 'Can only mask channel-independent quantities from radiance' )
      if ( ignore .and. .not. got(f_channels) ) call AnnounceError ( key, &
        & 'Specifying "ignore" without "channels" turns everything off!' )
    else
      ! Do some error checking for the optical depth issues
      if ( any(got((/ f_opticalDepth, f_opticalDepthCutoff /))) ) then
        if ( .not. all(got((/ f_opticalDepth, f_opticalDepthCutoff /))) ) &
          & call AnnounceError ( key, &
          & 'Must supply both opticalDepth and opticalDepthCutoff' )
        if ( qty%template%quantityType /= l_radiance .or. &
          &  opticalDepth%template%quantityType /= l_opticalDepth ) &
          & call AnnounceError ( key, 'Supplied quantity is not optical depth' )
        if ( qty%template%signal /= opticalDepth%template%signal .or. &
          &  qty%template%sideband /= opticalDepth%template%sideband ) &
          & call AnnounceError ( key, &
          & 'Optical depth does not match subsetted quantity' )
      end if
      ! Check for exactly one of height, ignore, instances, reset, surface
      if ( count ( got ( &
        & (/ f_height, f_ignore, f_instances, &
        & f_reset, f_resetAll, f_surface, f_where /) &
        & ) ) /= 1 ) &
          & call announceError ( key, &
            & 'Subset must be exactly one of height, ignore, instances, ' // &
            & 'surface, where or reset' )
    end if

    if ( vectors(mainVectorIndex)%globalUnit /= phyq_invalid ) then
      testUnit = vectors(mainVectorIndex)%globalUnit
    else
      testUnit = qty%template%unit
    end if
    if ( got ( f_minValue ) .and. ( minUnit /= testUnit ) ) &
      & call AnnounceError ( key, WrongUnits, f_minValue )
    if ( got ( f_maxValue ) .and. ( maxUnit /= testUnit ) ) &
      & call AnnounceError ( key, WrongUnits, f_maxValue )

    if ( resetAll ) then
      ! Simply Deallocatee the quantity's Mask
      call DestroyVectorQuantityMask ( qty, 'SetUpSubset' )
    else
      ! Pass all the field nodes to the workhorse subroutine
      ! Let it set the mask
      call ApplyMaskToQuantity( qty, rad, ptan, opticalDepth, a, &
        & opticalDepthCutoff, maxvalue, minValue, missingValue, &
        & heightRange, whereRange, &
        & ignore, reverse, additional, expandMask, reset, &
        & maskBit, heightNode, surfNode, instancesNode, channelsNode, got )
    endif

    if ( switchDetail(switches,'msk') > -1 ) then
      call output ( 'Dumping mask' )
      if ( qty%template%name /= 0 ) then
        call output ( ' for quantity ' )
        call display_string ( qty%template%name )
      end if
      call output ( '', advance='yes' )
      call DumpMask ( qty )
    end if
  end subroutine SetupSubset

  ! -------------------------------------------------- FlagCloud ---
  ! (Documented in wiki)
  ! Process flagCloud command. Set quantity mask if cloud contaminated.
  ! m_cloud is the default if no maskbit is given

  subroutine SetupFlagCloud ( key, vectors )
    
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_M, only: Expr, Getindexflagsfromlist
    use Init_Tables_Module, only: Field_First, Field_Last
    use Init_Tables_Module, only: F_Quantity, F_Ptanquantity, &
      & F_Height, F_Cloudheight, F_Channels, F_Cloudchannels, F_Cloudradiance, &
      & F_Cloudradiancecutoff, F_Mask
    use Init_Tables_Module, only: L_Radiance, L_Cloudinducedradiance
    use Intrinsic, only: Phyq_Dimensionless, Phyq_Pressure
    use MLSKinds, only: R8
    use Moretree, only: Get_Field_Id
    use VectorsModule, only: CreateMask, &
      & GetVectorQtyByTemplateIndex, SetMask, VectorValue_T, Vector_T, &
      & M_Cloud
    use Tree, only: Nsons, Subtree, Decoration

    integer, intent(in) :: KEY        ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    integer :: CHANNEL                ! Loop index
    integer :: CHANNELSNODE           ! Tree node for channels values
    integer :: FIELD                  ! Field type from tree
    integer :: GSON                   ! Tree node
    integer :: HEIGHT                 ! Loop counter
    integer :: HEIGHT1                 ! Loop counter
    integer :: CLOUDHEIGHTNODE             ! Tree node for cloudheight values
    integer :: HEIGHTNODE             ! Tree node for height values
    integer :: CLOUDHEIGHTUNIT             ! Unit for cloudheights command
    integer :: HEIGHTUNIT             ! Unit for heights command
    integer :: IND,IND1,J, I          ! Aarray indices
    integer :: INSTANCE               ! Loop counter
    integer :: MaskBit                ! Bits corresponding to Mask
    integer :: QUANTITYINDEX          ! Index
    integer :: SON                    ! Tree node
    integer :: STATUS                 ! Flag
    integer :: UNITS(2)               ! Units returned by expr
    integer :: VECTORINDEX            ! Index
    integer :: USETHISCHANNEL         ! cloud radiance channel

    real(r8), dimension(:), pointer :: THESEHEIGHTS ! Subset of heights
    real(r8) :: VALUE(2)              ! Value returned by expr
    real(r8) :: HeightVALUE(2)              ! Value returned by height expr
    real(r8) :: CHeightVALUE(2)              ! Value returned by cloudHeight expr
    real(r8) :: cloudRadianceCutOff              ! threshold for flagging cloud
    type (VectorValue_T), pointer :: QTY ! The quantity to mask
    type (VectorValue_T), pointer :: PTAN ! The ptan quantity if needed
    type (VectorValue_T), pointer :: cloudRadiance ! cloud radiances

    logical :: Got(field_first:field_last)   ! "Got this field already"
    logical, dimension(:), pointer :: CHANNELS ! Are we dealing with these channels
    logical :: DOTHISCLOUDHEIGHT           ! Flag
    logical :: DOTHISHEIGHT           ! Flag
    logical :: DOTHISCHANNEL          ! Flag
    logical :: ISCLOUD                ! Flag

    ! Executable code
    nullify ( channels, qty, ptan, cloudRadiance )
    got = .false.
    maskBit = m_cloud

    do j = 2, nsons(key) ! fields of the "Flagcloud" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        qty => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_ptanquantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        ptan => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_height )
        heightNode = son
      case ( f_cloudHeight )
        cloudHeightNode = son
      case ( f_channels )
        channelsNode = son
      case ( f_cloudChannels )
        if ( nsons(son) /= 2 ) call AnnounceError ( key, &
          & 'Too many cloudChannels' )
        call expr ( subtree(2,son), units, value )
        useThisChannel = nint ( value(1) )
      case ( f_cloudRadiance )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        cloudRadiance => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_cloudRadianceCutoff )
        call expr ( subtree (2, son), units, Value )
        cloudRadianceCutoff = value(1)
      case ( f_mask )
        maskBit = GetMaskBit ( son )
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! Quantity must be radiance
    if ( (.not. any(qty%template%quantityType==(/ l_radiance,l_cloudinducedradiance /))) &
      & .or. (cloudRadiance%template%quantityType /= l_cloudInducedRadiance &
      & .and. cloudRadiance%template%quantityType /= l_radiance)) &
      & call AnnounceError ( key, 'You supplied the wrong quantities to flagCloud' )

    if( got(f_channels)) then
      call Allocate_test ( Channels, qty%template%noChans, &
        & 'channels', ModuleName )
      call GetIndexFlagsFromList ( channelsNode, channels, status, &
        & lower=lbound(channels,1) )
    end if

    ! Preprocess the height stuff.
    heightUnit = phyq_dimensionless
    if ( got(f_height) ) then
      if (nsons(heightNode) .gt. 2) call AnnounceError ( key, &
        & 'Only one height range allowed', f_height )

      call expr ( subtree(2,heightNode), units, value )
      ! Make sure the range has non-dimensionless units -- the type
      ! checker only verifies that they're consistent.  We need to
      ! check each range separately, because the units determine the
      ! scaling of the values.
      if ( all(units == phyq_dimensionless) ) call announceError ( &
        & subtree(2,heightNode), WrongUnits, f_height )
      ! Check consistency of units -- all the same, or dimensionless. The
      ! type checker verifies the consistency of units of ranges, but not
      ! of array elements.
      do i = 1, 2
        if ( heightUnit == phyq_dimensionless ) then
          heightUnit = units(i)
        else if ( units(i) /= phyq_dimensionless .and. &
          &       units(i) /= heightUnit ) then
          call announceError ( heightNode, WrongUnits, f_height )
        end if
      end do

      ! convert pressure unit
      Heightvalue = value
      if ( heightUnit == phyq_pressure ) Heightvalue = -log10(value)

    end if
    ! Preprocess the height stuff.
    cloudheightUnit = phyq_dimensionless
    if ( got(f_cloudheight) ) then
      if (nsons(cloudheightNode) .gt. 2) call AnnounceError ( key, &
        & 'Can only supply one height range' ) 

      call expr ( subtree(2,cloudheightNode), units, value )
      if ( all(units == phyq_dimensionless) ) call announceError ( &
        & subtree(2,cloudheightNode), WrongUnits, f_cloudheight )
      do i = 1, 2
        if ( cloudheightUnit == phyq_dimensionless ) then
          cloudheightUnit = units(i)
        else if ( units(i) /= phyq_dimensionless .and. &
          &       units(i) /= cloudheightUnit ) then
          call announceError ( key, 'Problem with units', f_cloudheight )
        end if
      end do

      ! convert pressure unit
      CHeightvalue = value
      if ( cloudheightUnit == phyq_pressure ) CHeightvalue = -log10(value)          
    end if

    ! ----- finish checking ------

    ! Create the mask if it doesn't exist
    if ( .not. associated( qty%mask ) ) call CreateMask ( qty )

    ! Now loop over the instances
    do instance = 1, qty%template%noInstances
      theseHeights => ptan%values(:,instance)

      isCloud = .false.

      ! use the given height range to find cloud flag
      if( got(f_cloudheight) ) then 
        do height1 = 1,cloudRadiance%template%noSurfs
          doThisCloudHeight = .true.
          doThisCloudHeight = doThisCloudHeight .and. &
            & theseHeights(height1) > Cheightvalue(1)
          doThisCloudHeight = doThisCloudHeight .and. &
            & theseHeights(height1) < Cheightvalue(2)
          ind1 = useThisChannel + cloudRadiance%template%noChans*(height1-1)
          if (cloudRadianceCutoff > 0._r8) then
             if ( doThisCloudHeight .and. &
               & cloudRadiance%values ( ind1, instance ) > cloudRadianceCutoff) &
               & isCloud = .true.
          else
             if ( doThisCloudHeight .and. &
               & cloudRadiance%values ( ind1, instance ) < cloudRadianceCutoff) &
               & isCloud = .true.
          endif
        enddo
      endif
      do height = 1, qty%template%noSurfs
        doThisHeight = .true.
        doThisHeight = doThisHeight .and. &
          & theseHeights(height) > heightvalue(1)
        doThisHeight = doThisHeight .and. &
          & theseHeights(height) < heightvalue(2)

        ! determine cloud flag
        ! if cloud flag is not from a height range then do comparison at each minor frame
        if ( .not. got(f_cloudheight) ) then 
          ! use the same mif radiance to find cloud flag
          isCloud = .false.
          ind1 = useThisChannel + cloudRadiance%template%noChans*(height-1)
          if (cloudRadianceCutoff > 0._r8) then
            if ( cloudRadiance%values ( ind1, instance ) > cloudRadianceCutoff) &
            & isCloud = .true.
          else
            if ( cloudRadiance%values ( ind1, instance ) < cloudRadianceCutoff) &
            & isCloud = .true.
          endif
        endif

        ! assign cloud flag
        do channel = 1, qty%template%noChans
          doThisChannel = .true.
          if ( associated(channels) ) doThisChannel = channels(channel)
          if ( doThisChannel ) then
            ind = channel + qty%template%noChans*(height-1)
            if ( doThisHeight .and. isCloud )  &
              &     call SetMask ( qty%mask(:,instance), (/ ind /), &
              &     what=maskBit )
          end if                        ! do this channel
        end do                          ! Channel loop
      end do                            ! height loop
    end do                              ! Instance loop

    ! Tidy up
    if ( associated(channels) ) &
      & call Deallocate_test ( channels, 'channels', ModuleName )

  end subroutine SetupFlagCloud

  ! ---------------------------------------------- UpdateMask -----
  ! This is probably a nice operation
  ! Unfortunately, it lacks any documentation, so who knows
  ! what it does or whether it works properly?
  
  ! So be warned--
  ! d o c u m e n t   y o u r   c o d e
  subroutine UpdateMask ( key, vectors )

    use Init_Tables_Module, only: Field_First, Field_Last
    use Init_Tables_Module, only: F_Quantity, F_SourceQuantity, &
      & F_Operation, F_SourceMask, F_Mask
    use Init_Tables_Module, only: L_Invert, L_Copy, L_AndMasks, L_OrMasks
    use Moretree, only: Get_Field_Id
    use Tree, only: Nsons, Subtree, Decoration
    use VectorsModule, only: GetVectorQtyByTemplateIndex, VectorValue_T, &
      & Vector_T, CreateMask

    integer, intent(in) :: KEY          ! Tree node
    type (Vector_T), dimension(:) :: VECTORS

    ! Local variables
    logical :: GOT(field_first:field_last)   ! "Got this field already"
    integer :: SON                      ! Tree node
    integer :: GSON                     ! Another tree node
    integer :: MASKBIT                  ! Which bits
    integer :: SOURCEMASKBIT            ! Which bits from source
    integer :: J                        ! Loop counter
    integer :: FIELD                    ! A field id
    integer :: VECTORINDEX              ! In database
    integer :: QUANTITYINDEX            ! In database not vector
    integer :: OPERATION                ! What operation to perform
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: SOURCEQUANTITY

    ! Executable code
    ! First work out what we've been asked to do

    ! Executable code
    nullify ( sourceQuantity )
    got = .false.
    do j = 2, nsons(key) ! fields of the "subset" specification
      son = subtree(j, key)
      field = get_field_id(son)   ! tree_checker prevents duplicates
      if (nsons(son) > 1 ) gson = subtree(2,son) ! Gson is value
      select case ( field )
      case ( f_quantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex), quantityIndex)
      case ( f_sourceQuantity )
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        sourceQuantity => GetVectorQtyByTemplateIndex(vectors(vectorIndeX), quantityIndex)
      case ( f_mask )
        if ( .not. got(f_mask) ) maskBit = 0 ! clear default first time
        maskBit = GetMaskBit ( son, single=.true. )
      case ( f_sourceMask )
        if ( .not. got(f_sourceMask) ) sourceMaskBit = 0 ! clear default first time
        sourceMaskBit = GetMaskBit ( son, single=.true. )
      case ( f_operation )
        operation = decoration(gson)
      case default
        ! Shouldn't get here if the type checker worked
      end select
      got(field) = .true.
    end do ! j = 2, nsons(key)

    ! Special checking for invert case
    if ( operation == l_invert ) then
      if ( associated(sourceQuantity) ) call AnnounceError ( key, &
        & 'sourceQuantity not appropriate for invert updateMask' )
      if ( got(f_sourceMask) ) call AnnounceError ( key, &
        & 'sourceMask not appropriate for invert updateMask' )
    end if

    ! If we've got a source quantity then check it out.
    if ( associated ( sourceQuantity ) ) then
      if ( sourceQuantity%template%name /= quantity%template%name ) &
        & call AnnounceError ( key, 'sourceQuantity and quantity not compatible' )
      if ( .not. got(f_sourceMask) ) sourceMaskBit = maskBit
    else
      if ( operation /= l_invert ) then
        if ( .not. got(f_sourceMask) ) call AnnounceError ( key, &
          & 'sourceMask required for this operation' )
        if ( sourceMaskBit == maskBit ) call AnnounceError ( key, &
          & 'mask and sourceMask are the same information' )
      end if
      sourceQuantity => quantity
    end if

    ! If the result quantity has no mask then create one.
    if ( .not. associated ( quantity%mask ) ) call CreateMask ( quantity )
    ! Create one in the source too if it doesn't have one, potentially time wasting
    ! later on, but unlikely
    if ( .not. associated ( sourceQuantity%mask ) ) call CreateMask ( sourceQuantity )

    ! Now work out what to do.
    select case ( operation )
    case ( l_invert )
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       iand ( not ( iand ( ichar(quantity%mask), maskBit ) ), maskBit ) ) )
    case ( l_copy )
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( maskBit, 0, &
        &              iand ( ichar(sourceQuantity%mask), sourceMaskBit ) /= 0 ) ) )
    case ( l_andMasks )
      ! Note, because we 'invert' the logic in storing the mask
      ! this operation may look wrong at first glance
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( 0, maskBit, &
        &               iand ( ichar(sourceQuantity%mask), sourceMaskBit ) == 0 .and. &
        &               iand ( ichar(quantity%mask), maskBit ) == 0  ) ) )
    case ( l_orMasks )
      ! Note, because we 'invert' the logic in storing the mask
      ! this operation may look wrong at first glance
      quantity%mask = char ( &
        & ior ( iand ( ichar(quantity%mask), not(maskBit) ), & ! Other bits
        &       merge ( 0, maskBit, &
        &               iand ( ichar(sourceQuantity%mask), sourceMaskBit ) == 0 .or. &
        &               iand ( ichar(quantity%mask), maskBit ) == 0  ) ) )
    end select

  end subroutine UpdateMask
 
  ! ===================================== Private procedures ============

  ! --------------------------------------- AnnounceError -----------

  subroutine AnnounceError ( NODE, STRING, FIELDINDEX, ANOTHERFIELDINDEX )
    
    use Intrinsic, only: Field_Indices
    use Lexer_Core, only: Print_Source
    use Tree, only: Where
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_M, only: Output
    use String_Table, only: Display_String
    
    integer, intent(in) :: NODE
    character(len=*), intent(in) :: STRING
    integer, intent(in), optional :: FIELDINDEX, ANOTHERFIELDINDEX ! f_...
    
    call output ( '***** At ' )
    call print_source ( where(node) )
    call output ( ', SubsetModule complained: ', advance='yes' )
    call output ( trim(string), advance='yes' )
    if ( present ( fieldIndex ) ) then
      call output ( 'Offending field(s): ' )
      call display_string ( field_indices(fieldIndex) )
      if ( present ( anotherFieldIndex ) ) then
        call output ( ', ' )
        call output ( field_indices(anotherFieldIndex) )
      end if
      call output ( '', advance='yes' )
    end if
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Problem in SubsetModule, see above' )
  end subroutine AnnounceError

  ! -------------------------------------- GetMaskBit -------------

  function GetMaskBit ( node, single ) result ( maskBit )
    use Tree, only: Nsons, Decoration, Subtree
    use Init_Tables_Module, only: L_Cloud, L_Fill, L_Full_Derivatives, &
      & L_Ignore, L_Linalg, L_Spare, L_Tikhonov
    use VectorsModule, only: M_Cloud, M_Fill, M_Fullderivatives, M_Linalg, &
      & M_Ignore, M_Spare, M_Tikhonov

    ! This routine parses the mask field in an l2cf command
    ! Dummy arguments
    integer, intent(in) :: NODE         ! Tree node
    logical, optional, intent(in) :: SINGLE ! Raise error if more than one
    integer :: MASKBIT                  ! Result
    integer :: MASK                     ! Tree decoration
    integer :: I
    if ( present ( single ) ) then
      if ( single .and. nsons(node) > 2 ) &
        & call AnnounceError ( node, 'Only one mask expected' )
    end if
    maskBit = 0
    do i = 2, nsons(node)
      mask = decoration(subtree(i,node))
      select case ( mask )
      case ( l_Cloud )
        maskBit = ior(maskBit, m_Cloud)
      case ( l_fill )
        maskBit = ior(maskBit, m_fill)
      case ( l_full_derivatives )
        maskBit = ior(maskBit, m_fullDerivatives)
      case ( l_ignore )
        maskBit = ior(maskBit, m_Ignore)
      case ( l_linAlg )
        maskBit = ior(maskBit, m_linAlg)
      case ( l_Spare )
        maskBit = ior(maskBit, m_Spare)
      case ( l_Tikhonov )
        maskBit = ior(maskBit, m_Tikhonov)
      end select
    end do
  end function GetMaskBit

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SubsetModule
 
! $Log$
! Revision 2.39  2019/10/16 20:55:04  pwagner
! Subset command may take a MissingValue field
!
! Revision 2.38  2017/07/10 18:50:00  pwagner
! Transfer may /expandMask to all masking bits; may /skipValues to transfer only mask; Fill may replaceMissingValue=
!
! Revision 2.37  2017/02/22 01:23:36  pwagner
! New /resetAll switch clears all Masking bits
!
! Revision 2.36  2016/05/19 23:17:56  pwagner
! Corrected what looked like an error in indexing
!
! Revision 2.35  2016/05/04 17:52:00  pwagner
! Tried harder to make the hightRange mechanism work as advertised
!
! Revision 2.34  2014/09/05 01:23:40  vsnyder
! Remove declaration of unused parameter and USE name
!
! Revision 2.33  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.32  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.31  2014/02/28 01:15:27  vsnyder
! Remove TYPE argument from calls to EXPR because the value wasn't used
!
! Revision 2.30  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.29  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.28  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.27  2013/01/18 01:45:57  pwagner
! ApplyMaskToQuantity can set mask based on values of 'a' satisfying condition 'where'
!
! Revision 2.26  2012/10/22 18:14:11  pwagner
! Many Subset operations now available in Fill
!
! Revision 2.25  2012/05/03 19:50:11  vsnyder
! More error checking for mask from radiance
!
! Revision 2.24  2012/05/01 22:22:26  vsnyder
! Mask minor-frame quantity by reference to radiance quantity.  Move some
! dump generation to TrueList in MLSStrings.
!
! Revision 2.23  2012/03/28 00:54:41  vsnyder
! Better error message for wrong units
!
! Revision 2.22  2011/05/09 18:26:45  pwagner
! Converted to using switchDetail
!
! Revision 2.21  2011/03/15 22:54:18  pwagner
! /reverse flag reverses which elements are masked by Subset command
!
! Revision 2.20  2010/04/28 16:24:43  pwagner
! May specify instances range in Subset
!
! Revision 2.19  2010/04/22 23:38:11  pwagner
! Added new Ignore masking bit
!
! Revision 2.18  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.17  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.16  2006/04/04 15:59:48  pwagner
! NAG caught us mixing variable types in min, max
!
! Revision 2.15  2006/04/03 23:54:33  livesey
! Added range guards
!
! Revision 2.14  2006/04/03 23:51:42  livesey
! Added f_surface capability to subset
!
! Revision 2.13  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.12  2004/05/28 00:57:50  vsnyder
! Move GetIndexFlagsFromList from MoreTree to Expr_m
!
! Revision 2.11  2003/08/06 17:24:16  livesey
! Essentially cosmetic change
!
! Revision 2.10  2003/05/14 22:02:08  dwu
! fix a bug in flagcloud
!
! Revision 2.9  2003/05/12 18:55:01  dwu
! allow cloud flag for cloudinducedradiance
!
! Revision 2.8  2003/04/14 22:21:26  dwu
! make cloudradiancecutoff sign dependent: it is 'greater than' if positive; 'less than' if negative
!
! Revision 2.7  2003/04/10 18:30:57  dwu
! make m_cloud as default in FlagCloud
!
! Revision 2.6  2003/04/08 23:12:18  dwu
! add m_cloud in setupFlagCloud
!
! Revision 2.5  2003/04/04 22:01:46  livesey
! Added UpdateMask
!
! Revision 2.4  2003/03/19 01:58:36  livesey
! Added nullify on signals
!
! Revision 2.3  2003/03/07 03:17:12  livesey
! Added RestrictRange
!
! Revision 2.2  2003/03/06 00:47:04  livesey
! Added logging information
!
