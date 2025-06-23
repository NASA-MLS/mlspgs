! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================

module L2PCBins_m

  implicit none
  private
  public :: FindMatchForL2PCQ, FlushLockedBins, SelectL2PCBins

  ! This array is used to keep track of which bins to use for each (side)band.
  integer, dimension(:,:), pointer, private, save :: lockedBins => NULL()

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ----------------------------------------  FindMatchForL2PCQ  -----
  subroutine FindMatchForL2PCQ ( l2pcQ, fmConf, FwdModelIn, FwdModelExtra, &
    & stateQ, foundInFirst )

    use Forwardmodelconfig, only: Forwardmodelconfig_T
    use ForwardmodelVectortools, only: Getquantityforforwardmodel
    use Intrinsic, only: L_Fieldazimuth, L_Fieldelevation, L_Fieldstrength, &
      & L_Temperature, L_Tscat, L_Vmr
    use ManipulateVectorquantities, only: Dofgridsmatch, Dohgridsmatch, &
      & Dovgridsmatch
    use Molecules, only: Isextinction
    use Vectorsmodule, only: VectorValue_T, Vector_T

    type (VectorValue_T), intent(in) :: L2PCQ ! Quantity to search for
    type (ForwardModelConfig_T), intent(in) :: fmConf ! Forward model config
    type (Vector_T), intent(in), target :: FWDMODELIN ! State vector
    type (Vector_T), intent(in), target :: FWDMODELEXTRA ! Extra state vector
    type (VectorValue_T), pointer :: STATEQ ! Result
    logical, intent(out) :: FOUNDINFIRST ! If set, found in first vector
    ! This routine looks through the supplied state vectors and finds a match
    ! for the supplied quantity from the l2pc state vector.  It then goes on to
    ! test that the HGrids and VGrids for the two quantities match appropriately.

    ! Local variables
    integer :: QTY, VEC
    type (Vector_T), pointer :: V

    ! Executable code
    foundInFirst = .false.
    ! We're going to be fairly picky about this.  I don't want to get into the
    ! habit of blindly accepting quantities that might not be appropriate.
    ! Therefore, I have a list of 'acceptable' quantity types.
    select case ( l2pcQ%template%quantityType )
    case ( l_fieldStrength, l_fieldAzimuth, l_fieldElevation, l_temperature, &
      &    l_TScat )
      ! This is for quantities that are 'easy to get'
      stateQ => GetQuantityForForwardModel ( FwdModelIn, FwdModelExtra,&
        & quantityType = l2pcQ%template%quantityType, config=fmConf, &
        & foundInFirst = foundInFirst, noError=.true. )
    case ( l_vmr )
      ! Here we may need to be a little more intelligent
      if ( .not. isExtinction(l2pcQ%template%molecule) ) then
        stateQ => GetQuantityForForwardModel ( FwdModelIn, FwdModelExtra,&
          & quantityType = l_vmr, config=fmConf, &
          & molecule = l2pcQ%template%molecule, &
          & foundInFirst = foundInFirst, noError=.true., matchQty=l2pcQ )
      else
        searchLoop: do vec = 1, 2
          ! Point to appropriate vector
          if ( vec == 1 ) then
            v => FwdModelIn
          else
            v => FwdModelExtra
          end if
          ! Skip this vector if it doesn't exist
          if ( .not. associated ( v ) ) cycle searchLoop
          ! Loop over this vector
          do qty = 1, size ( v%quantities )
            stateQ => v%quantities(qty)
            if ( stateQ%template%quantityType == l_vmr .and. &
              &  isExtinction(l2pcQ%template%molecule) .and. &
              &  stateQ%template%radiometer == l2pcQ%template%radiometer ) then
              if ( DoFGridsMatch ( l2pcQ, stateQ ) ) exit searchLoop
            end if
          end do
        end do searchLoop
        foundInFirst = ( vec == 1 )
      end if
    case default
      nullify ( stateQ )
    end select

    ! Now check that these match.
    if ( associated ( stateQ ) ) then
      if ( .not. DoVGridsMatch ( stateQ, l2pcQ ) ) stateQ => NULL()
    end if
    if ( associated ( stateQ ) ) then
      if ( .not. DoHGridsMatch ( stateQ, l2pcQ, spacingOnly=.true. ) ) stateQ => NULL()
    end if

    if ( .not. associated ( stateQ ) ) foundInFirst = .false.

  end subroutine FindMatchForL2PCQ

  ! --------------------------------------------  FlushLockedBins  -----
  subroutine FlushLockedBins
    use Allocate_Deallocate, only: Deallocate_test
    use L2PC_m, only: FlushL2PCBins

    ! This could set them to zero again, but I think I'll deallocate here, just
    ! because otherwise I'd have to write a separate deallocate routine.
    call deallocate_test ( lockedBins, 'lockedBins', moduleName )
    call FlushL2PCBins
  end subroutine FlushLockedBins

  ! ---------------------------------------------  SelectL2PCBins  -----
  subroutine SelectL2PCBins ( fmConf, FwdModelIn, FwdModelExtra, &
    & radiance, sideband, maf, &
    & l2pcBins, sidebandStart, sidebandStop, sidebandStep )

    ! Selection criteria
    ! 1.  Signal in L2PCdatabase(bin) == signal
    ! 2.  L2PCdatabase(bin) matches some binSelectors(fmConf%binSelectors)
    ! 3.  For folded-sideband case, prefer folded-sideband bin
    ! 4.  From the eligible ones, pick the minimum cost, which depends
    !     upon selectorType

    use Allocate_Deallocate, only: Allocate_Test
    use Forwardmodelconfig, only: Forwardmodelconfig_T
    use Intrinsic, only: L_Fieldazimuth, L_Fieldelevation, L_Fieldstrength, &
      & L_Latitude, L_Namefragment, L_Sza, L_Temperature, L_Tscat, L_Vmr, L_Zeta
    use L2pc_M, only: Binselectors, Binselector_T, L2pcDatabase
    use ManipulateVectorquantities, only: Findoneclosestinstance
    use HyperSlabs, only: Essentiallyequal
    use MLSKinds, only: R8
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_M, only: Getsidebandloop, Getsignalname, Signals
    use MLSStringlists, only: Switchdetail
    use Output_M, only: Output
    use QuantityTemplates, only: QuantityTemplate_T
    use String_Table, only: Display_String, Index, Len
    use Toggles, only: Switches
    use Vectorsmodule, only: GetVectorquantitybytype, ValidateVectorquantity, &
      & VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(in) :: FMCONF
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type (VectorValue_T), intent(in) :: RADIANCE ! The radiance we're after
    integer, intent(in) :: SIDEBAND ! Which sideband (see below)
    integer, intent(in) :: MAF                      ! MAF index
    integer, dimension(-1:1), intent(out) :: L2PCBINS ! Result
    integer, intent(out), optional :: SIDEBANDSTART ! Resulting loop indices
    integer, intent(out), optional :: SIDEBANDSTOP
    integer, intent(out), optional :: SIDEBANDSTEP
    ! We supply the sideband separately because when the Linearized model is
    ! being invoked by the hybrid model, the sideband in radiance is not the
    ! one we really want.
    ! When invoked for TScat from the full forward model, it's already
    ! decided which sideband it's working on.

    ! Local variables
    character(132) :: SignalName
    integer :: BIN                    ! Loop counter
    integer :: L2PCINSTANCE           ! Instance index
    integer :: MAF1                   ! Subset limit
    integer :: MAFN                   ! Subset limit
    integer :: MYMAF                  ! A loop counter
    integer :: NOBINS                 ! Number of l2pc bins
    integer :: NOSELECTORS            ! A loop limit
    integer :: SELECTOR               ! Bin selector index
    integer :: SIGNAL                 ! Signal index
    integer :: STATEINSTANCE          ! Instance index
    integer :: tmpSideband            ! Loop counter
    logical :: SPLIT                  ! Need a split calculation
    logical :: FOUNDINFIRST           ! Flag, not used
    real(r8) :: THISCOST              ! Interim cost

    integer, dimension(1) :: S1, S2   ! Results of minloc
    logical, dimension(-1:1,size(l2pcDatabase)) :: POSSIBLE ! Flags for each bin
    real(r8), dimension(size(l2pcDatabase)) :: COST ! Cost of each bin
    real(r8), dimension(-1:1) :: BESTCOST
    type (QuantityTemplate_T), pointer :: BINRAD ! Quantity template
    type (BinSelector_T), pointer :: SEL         ! One bin selector
    type (VectorValue_T), pointer :: L2PCQ       ! A quantity in the l2pc
    type (VectorValue_T), pointer :: STATEQ      ! A state vector quantity

    ! Executable code

    ! Firstly, if we're in locked bins mode, then just return the locked bin
    signal = radiance%template%signal
    noBins = size ( l2pcDatabase )

    if ( signal == 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "No signal available when selecting L2PC bin" )

    ! Now special code for dealing with the locked bins case
    if ( fmConf%lockBins ) then
      if ( .not. associated ( lockedBins ) ) &
        & call allocate_test ( lockedBins, 1, size(signals), 'lockedBins', &
          & moduleName, low1=-1, fill=0 )
      if ( any ( lockedBins ( :, signal ) /= 0 ) ) then
        l2pcBins = lockedBins ( :, signal )
        ! If got folded, use that by preference (3rd argument below).
        ! Assert that sidebands are present if needed here.
        if ( present(sidebandStart) ) &
          & call GetSidebandLoop ( signal, sideband, &
            & (lockedBins(0,signal)==0), sidebandStart, sidebandStop, &
            & sidebandStep )
        return
      end if
    end if

    ! Code will only get to here if the bins are unlocked, or will be locked
    ! once we've decided on them.

    ! Setup the arrays we need, possible is a flag to indicate (for
    ! each sideband) whether a bin is even worth considering, cost is
    ! the associated cost.

    ! Get a first cut at the 'possible' flag
    possible = .false.
    do bin = 1, noBins
      binRad => l2pcDatabase(bin)%j%row%vec%quantities(1)%template
      possible ( binRad%sideband, bin ) = ( binRad%signal == signal )
    end do
    ! Set the cost to zero by default
    cost = 0.0_r8

    ! Setup some stuff we'll need from time to time
    if ( fmConf%lockBins ) then
      maf1 = 1 + radiance%template%noInstancesLowerOverlap
      mafN = radiance%template%noInstances - &
        & radiance%template%noInstancesUpperOverlap
    else
      maf1 = maf
      mafN = maf
    end if

    ! OK, now we have to loop over the bin selectors for each bin and
    ! apply the rules they give.
    noSelectors = size ( fmConf%binSelectors )
    do bin = 1, noBins
      if ( .not. any ( possible ( :, bin ) ) ) cycle
      binRad => l2pcDatabase(bin)%j%row%vec%quantities(1)%template
      do selector = 1, noSelectors
        sel => binSelectors ( fmConf%binSelectors(selector) )
        select case ( sel%selectorType )
        case ( l_nameFragment )
          if ( l2pcDatabase(bin)%name < 1 ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'l2pc db name missing from string table' )
          if ( sel%nameFragment < 1 ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'l2pc selector name fragment missing from string table' )
          if ( len(sel%nameFragment) > 0 ) &
              & possible ( :, bin ) = possible ( :, bin ) .and. &                                 
                & index ( l2pcDatabase(bin)%name, sel%nameFragment, &
                  & caseless=.true., strip=.true. ) /= 0
        case ( l_latitude )
          ! When we say latitude, we really mean an empirical phi.
          thisCost = sqrt ( sum ( &
            & (NormalizePhi(radiance%template%phi(1,maf1:mafN)) - &  
            &  NormalizePhi(binRad%phi(1,1)) )**2 ) / (mafN-maf1+1) ) / &
            & sel%cost
          if ( sel%exact ) then
            possible ( :, bin ) = possible ( :, bin ) .and. &
              & EssentiallyEqual ( thisCost, 0.0_r8 )
          else
            cost(bin) = cost(bin) + thisCost
          end if
        case ( l_sza )
          thisCost = sqrt ( sum ( &
            & ( radiance%template%solarZenith(1,maf1:mafN) - &
            &   binRad%solarZenith(1,1) )**2 ) / &
            & (mafN-maf1+1) ) / sel%cost
          if ( sel%exact ) then
            possible ( :, bin ) = possible ( :, bin ) .and. &
              & EssentiallyEqual ( thisCost, 0.0_r8 )
          else
            cost(bin) = cost(bin) + thisCost
          end if
        case ( l_temperature, l_fieldAzimuth, l_fieldElevation, l_fieldStrength, &
          &    l_TScat, l_vmr )
          ! This one involves matching elements of xStar with x.
          if ( sel%selectorType == l_vmr ) then
            l2pcQ => GetVectorQuantityByType ( &
              & l2pcDatabase(bin)%j%col%vec, quantityType=l_vmr, &
              & molecule=sel%molecule, noError=.true. )
          else
            l2pcQ => GetVectorQuantityByType ( &
              & l2pcDatabase(bin)%j%col%vec, quantityType=sel%selectorType, &
              & noError=.true. )
          end if
          if ( .not. associated ( l2pcQ ) ) then
            possible ( :, bin ) = .false.
          else
            if ( .not. ValidateVectorQuantity ( l2pcQ,&
              & verticalCoordinate = (/ l_zeta /) ) ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Expected zeta coordinate for quantity in binSelector' )
            ! Find the relevant corresponding quantity in the state vector
            call FindMatchForL2PCQ ( l2pcQ, fmConf, fwdModelIn, &
              & fwdModelExtra, stateQ, foundInFirst )
            ! If we've got both of them make sure they match
            if ( associated(stateQ) .and. associated(l2pcQ) ) then
              ! OK, identify height range
              if ( all ( binSelectors(selector)%heightRange > 0.0 ) ) then
                s1 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                  & log10(binSelectors(selector)%heightRange(1) ) ) )
                s2 = minloc ( abs ( stateQ%template%surfs(:,1) + &
                  & log10(binSelectors(selector)%heightRange(2) ) ) )
              else
                s1 = 1
                s2 = stateQ%template%noSurfs
              end if
              ! Here we'll just compare the central profile in the
              ! l2pc with the state profile closest to each maf.
              l2pcInstance = l2pcQ%template%noInstances/2 + 1
              thisCost = 0.0_r8
              do myMaf = maf1, mafN
                ! Only compare the closest profile to the maf
                stateInstance = FindOneClosestInstance ( stateQ, radiance, myMaf )
                thisCost = thisCost + sum ( &
                  & ( stateQ%values ( s1(1):s2(1), stateInstance ) - &
                  &   l2pcQ%values  ( s1(1):s2(1), l2pcInstance  ) ) **2 )
              end do
              ! If we require an exact match set the possible flag,
              ! otherwise just report our cost.
              if ( sel%exact ) then
                possible ( :, bin ) = possible ( :, bin ) .and. &
                  & EssentiallyEqual ( thisCost, 0.0_r8 )
              else
                cost ( bin ) = cost ( bin ) + sqrt ( thisCost / &
                  &  ( ( s2(1)-s1(1)+1 ) * ( mafN - maf1 + 1 ) ) ) / sel%cost
              end if
            end if
          end if
        end select                  ! Bin selector type
      end do                        ! Loop over selectors
    end do                          ! Loop over bins

    if ( switchDetail ( switches, 'binsel' ) > -1 ) then
      call output ( 'Choosing bin for ' )
      call GetSignalName ( signal, signalName, sideband=sideband )
      call output ( trim(signalName), advance='yes' )
      do bin = 1, noBins
        if ( any ( possible ( :, bin ) ) ) then
          call output  ( 'Candidate: ' )
          call display_string ( l2pcDatabase(bin)%name, strip=.true. )
          call output ( ' cost = ' )
          call output ( cost(bin), advance='yes' )
        end if
      end do
    end if

    ! OK, now we've surveyed the bins, let's cut things down.
    ! When computing folded radiances I'll always choose folded bins
    ! over unfolded ones, even if they cost more.  I think it's
    ! unlikely that this will really be an issue anyway.
    split = .false.
    if ( sideband == 0 ) then
      ! If we've got a match for the folded case, forget the others.
      if ( any ( possible(0,:) ) ) then
        where ( possible(0,:) ) 
          possible(-1,:) = .false.
          possible(1,:) = .false.
        end where
      else
        split = .true.
      end if
    end if

    if ( present(sidebandStart) ) then
      ! Work out the range of the sideband loop
      call GetSidebandLoop ( signal, sideband, split, &
        & sidebandStart, sidebandStop, sidebandStep )

      ! Choose the bin(s)
      bestCost = huge ( cost(1) )
      l2pcBins = 0
      do bin = 1, noBins
        do tmpSideband = sidebandStart, sidebandStop
          if ( possible(tmpSideband,bin) .and. &
            & cost(bin) < bestCost(tmpSideband) ) then
            bestCost(tmpSideband) = cost(bin)
            l2pcBins(tmpSideband) = bin
          end if
        end do
      end do

      ! Check that we've got the bins we need
      if (any(l2pcBins ((/sidebandStart,sidebandStop/)) == 0 )) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to find l2pc bins to match request' )

    else

      ! Choose the bin(s)
      bestCost = huge ( cost(1) )
      l2pcBins = 0
      do bin = 1, noBins
        if ( possible(sideband,bin) .and. &
          & cost(bin) < bestCost(sideband) ) then
          bestCost(sideband) = cost(bin)
          l2pcBins(sideband) = bin
        end if
      end do

    end if

    ! Record this in the 'locked bins' information
    if ( fmConf%lockBins ) lockedBins ( :, signal ) = l2pcBins

  contains

    ! ---------------------------------------------  NormalizePhi  -----
    elemental real (r8) function NormalizePhi ( phi ) result ( lat )
      ! This does an approximate phi to latitude conversion
      real (r8), intent(in) :: PHI      ! Input geodetic angle
      ! Executable code
      lat = modulo ( phi, 360.0_r8 )
      if ( (lat > 90.0) .and. ( lat <= 270.0 ) ) then
        lat = 180.0 - lat
      else if ( lat > 270.0 ) then
        lat = lat - 360.0
    ! else if ( lat <= 90/0 ) then
      ! nothing need be done
      end if
    end function NormalizePhi

  end subroutine SelectL2PCBins

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2PCBins_m

! $Log$
! Revision 2.10  2017/11/03 20:58:32  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.9  2011/11/11 00:42:06  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.8  2011/06/16 20:21:00  vsnyder
! Require a signal in the Radiance argument to SelectL2PCBins
!
! Revision 2.7  2011/05/09 17:48:37  pwagner
! Converted to using switchDetail
!
! Revision 2.6  2010/08/27 23:42:10  vsnyder
! Make SidebandStart etc. optional
!
! Revision 2.5  2010/05/19 17:52:53  pwagner
! Removed unused stuff
!
! Revision 2.4  2010/05/07 02:24:47  vsnyder
! Add strip=.true to a reference to string_table%index
!
! Revision 2.3  2010/04/30 23:57:07  vsnyder
! Add TScat.  Remove StateQ from SelectL2PCbins calling sequence.  Use
! INDEX from String_Table instead of fetching strings for intrinsic INDEX.
!
! Revision 2.2  2010/04/17 01:45:16  vsnyder
! Simplify some stuff
!
! Revision 2.1  2010/04/13 01:42:01  vsnyder
! Initial commit to move stuff here from LinearizedForwardModel_m
!
