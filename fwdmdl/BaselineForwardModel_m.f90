! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module BaselineForwardModel_m

  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, FORWARDMODELSTATUS_T
  use MLSCommon, only: RP
  use MLSSignals_m, only: SIGNALS, SIGNAL_T
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQUANTITYBYTYPE, &
    & VALIDATEVECTORQUANTITY
  use MatrixModule_1, only: MATRIX_T, FINDBLOCK, CREATEBLOCK
  use MatrixModule_0, only: SPARSIFY, MATRIXELEMENT_T, M_ABSENT, M_BANDED
  use Intrinsic, only: L_BASELINE, L_PTAN, L_NONE, L_RADIANCE, L_INTERMEDIATEFREQUENCY
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, &
    & MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use MLSNumerics, only: HUNT
  use Dump_0, only: DUMP

  ! This module contains a special forward model for baseline related effects.

  implicit none
  private
  public :: BaselineForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! ======================================== BaselineForwardModel ======

  subroutine BaselineForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    & FwdModelOut, oldIFM, fmStat, jacobian )

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Local parameters
    character(len=*), parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

    ! Local variables

    logical :: BSLINFIRST               ! Set if baseline in FwdModelIn
    logical :: PTANINFIRST              ! Set if ptan in FwdModelIn

    integer :: SIGINDEX                 ! Index into fwdModelConf%signals
    integer :: MAF                      ! Major frame index
    integer :: INSTANCE                 ! Instance index
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOCHANS                  ! Number of channels in radiance
    integer :: NOBSLCHANS               ! Number of channels in baseline
    integer :: ROW                      ! Element of jacobian
    integer :: ROWBLOCK                 ! Location in jacobian
    integer :: COLBLOCK                 ! Location in jacobian
    integer :: MIF                      ! Loop counter
    integer :: CHAN                     ! Loop counter
    integer :: INSTLOW                  ! Array limit
    integer :: INSTHI                   ! Array limit
    integer :: MM1                      ! MIF-1
    integer :: STATUS                   ! From allocates etc.

    real(rp) :: RAD                     ! One radiance
    real(rp) :: GRADIENT                ! A gradient

    integer, dimension(:), pointer :: CHAN0 ! Index into bsl freqs.
    integer, dimension(:), pointer :: CHAN1 ! Index into bsl freqs.
    integer, dimension(:), pointer :: INST0 ! Instance index into baseline
    integer, dimension(:), pointer :: INST1 ! Instance index into baseline
    integer, dimension(:), pointer :: SURF0 ! Surface index into baseline
    integer, dimension(:), pointer :: SURF1 ! Surface index into baseline
    integer, dimension(:), pointer :: SURF0m ! surf0 - 1
    integer, dimension(:), pointer :: SURF1m ! surf1 - 1

    real (rp), dimension(:), pointer :: CHANWT0 ! Weight for lower point
    real (rp), dimension(:), pointer :: CHANWT1 ! Weight for upper point
    real (rp), dimension(:), pointer :: INSTWT0 ! Weight for lower point
    real (rp), dimension(:), pointer :: INSTWT1 ! Weight for upper point
    real (rp), dimension(:), pointer :: SURFWT0 ! Weight for lower point
    real (rp), dimension(:), pointer :: SURFWT1 ! Weight for upper point
    real (rp), dimension(:), pointer :: SURFWT0PRIME ! d[SurfWt0]/d[ptan]
    real (rp), dimension(:), pointer :: SURFWT1PRIME ! d[SurfWt1]/d[ptan]
    real (rp), dimension(:,:), pointer :: KBIT2 ! Part of derivatives
    real (rp), dimension(:,:,:), pointer :: KBIT ! Part of derivatives

    type (VectorValue_T), pointer :: RADIANCE ! The radiance quantity
    type (VectorValue_T), pointer :: BASELINE ! The baseline quantity
    type (VectorValue_T), pointer :: PTAN ! The tangent pressure quantity

    type (MatrixElement_T), pointer :: JBLOCK     ! A block from the jacobian

    type (Signal_T), pointer :: SIGNAL ! The current signal being handled

    ! Executable code -------------------------------------------------------
    
    if (.not. fwdModelConf%do_Baseline ) return
    nullify ( chan0, chan1, inst0, inst1, surf0, surf1 )
    nullify ( chanWt0, chanWt1, instWt0, instWt1, surfWt0, surfWt1 )
    nullify ( surf0m, surf1m, surfWt0Prime, surfWt1Prime )
    nullify ( kBit )

    maf = fmStat%maf
    do sigIndex = 1, size(fwdModelConf%signals)
      signal => fwdModelConf%signals(sigIndex)

      ! Look for the radiance quantity for this signal
      radiance => GetVectorQuantityByType ( fwdModelOut, quantityType=l_radiance, &
        & signal=signal%index, sideband=signal%sideband )
      
      ! Set some dimensions
      noChans = radiance%template%noChans
      noMIFs = radiance%template%noSurfs

      ! Now work out which baseline quantity we want
      ! First look for one for this band
      baseline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_baseline, signal=signal%index, sideband=signal%sideband,&
        & noError=.true., foundInFirst=bslInFirst )

      ! If we can't find one, look for one for this radiometer instead,
      ! if that fails, raise an error
      if ( .not. associated(baseline) ) &
        & baseline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_baseline, radiometer=signal%radiometer, foundInFirst=bslInFirst )

      noBslChans = baseline%template%noChans

      ! Get ptans, we'll need these for interpolation
      ptan => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra, &
        & quantityType = l_ptan, &
        & instrumentModule = radiance%template%instrumentModule,&
        & foundInFirst=ptanInFirst )

      if (present(jacobian) .and. (bslInFirst .or. ptanInFirst) ) then
        rowBlock = FindBlock ( jacobian%row, radiance%index, maf )
        fmStat%rows(rowBlock) = .true.
      endif

      ! Now check the validity of the quantities we've been given
      if ( .not. ValidateVectorQuantity(baseline, stacked=.true., coherent=.true., &
        & regular=.true., &
        & frequencyCoordinate=(/ l_none, l_intermediateFrequency/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & InvalidQuantity//'baseline' )
      if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
        & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'ptan' )

      ! Now we identify the elements of baseline that will affect each radiance
      ! and compute weights

      ! Horizontal coordinate --------------------
      call Allocate_test ( inst0, noMIFs, 'inst0', ModuleName )
      call Allocate_test ( inst1, noMIFs, 'inst1', ModuleName )
      call Allocate_test ( instWt0, noMIFs, 'instWt0', ModuleName )
      call Allocate_test ( instWt1, noMIFs, 'instWt1', ModuleName )
      
      call Hunt ( baseline%template%phi(1,:), ptan%template%phi(:,maf), inst0 )
      inst1 = min ( inst0+1, baseline%template%noInstances )
      where ( inst1 /= inst0 )
        instWt1 = ( ptan%template%phi(:,maf) - baseline%template%phi(1,inst0) ) / &
          & ( baseline%template%phi(1,inst1) - baseline%template%phi(1,inst0) ) 
      elsewhere
        instWt1 = 0.0
      end where
      instWt0 = 1 - instWt1
      instWt1 = max(min(instWt1,1.0_rp),0.0_rp)
      instWt0 = max(min(instWt0,1.0_rp),0.0_rp)

      ! Vertical coordinate ---------------------
      call Allocate_test ( surf0, noMIFs, 'surf0', ModuleName )
      call Allocate_test ( surf1, noMIFs, 'surf1', ModuleName )
      call Allocate_test ( surf0m, noMIFs, 'surf0m', ModuleName )
      call Allocate_test ( surf1m, noMIFs, 'surf1m', ModuleName )
      call Allocate_test ( surfWt0, noMIFs, 'surfWt0', ModuleName )
      call Allocate_test ( surfWt1, noMIFs, 'surfWt1', ModuleName )
      
      call Hunt ( baseline%template%surfs(:,1), ptan%values(:,maf), surf0 )
      surf1 = min ( surf0+1, baseline%template%noSurfs )
      where ( surf1 /= surf0 )
        surfWt1 = ( ptan%values(:,maf) - baseline%template%surfs(surf0,1) ) / &
          & ( baseline%template%surfs(surf1,1) - baseline%template%surfs(surf0,1) ) 
      elsewhere
        surfWt1 = 0.0
      end where
      surfWt0 = 1 - surfWt1
      surfWt1 = max(min(surfWt1,1.0_rp),0.0_rp)
      surfWt0 = max(min(surfWt0,1.0_rp),0.0_rp)
      surf0m = surf0 - 1
      surf1m = surf1 - 1

      ! Frequency coordinate -------------------
      call Allocate_test ( chan0, noChans, 'chan0', ModuleName )
      call Allocate_test ( chan1, noChans, 'chan1', ModuleName )
      call Allocate_test ( chanWt0, noChans, 'chanWt0', ModuleName )
      call Allocate_test ( chanWt1, noChans, 'chanWt1', ModuleName )
      
      if ( associated ( baseline%template%frequencies ) ) then
        call Hunt ( baseline%template%frequencies, &
          & signal%frequencies+signal%centerFrequency, chan0 )
      else
        chan0 = 1
      end if
      chan1 = min ( chan0+1, noBslChans )
      if ( associated ( baseline%template%frequencies ) ) then
        where ( chan1 /= chan0 )
          chanWt1 = ( signal%frequencies+signal%centerFrequency - &
            & baseline%template%frequencies(chan0) ) / &
            & ( baseline%template%frequencies(chan1) - &
            &   baseline%template%frequencies(chan0) ) 
        elsewhere
          chanWt1 = 0.0
        end where
      else
        chanWt1 = 0.0
      endif
      chanWt0 = 1 - chanWt1
      chanWt1 = max(min(chanWt1,1.0_rp),0.0_rp)
      chanWt0 = max(min(chanWt0,1.0_rp),0.0_rp)

      ! -------------------------------------------------------------------
      ! Do this in a loop to avoid large array temps
      do mif = 1, noMIFs
        do chan = 1, noChans
          rad = radiance%values ( chan + noChans*(mif-1), maf )
          rad = rad + chanWt0(chan) * surfWt0(mif) * instWt0(mif) * &
            & baseline%values(chan0(chan)+noBslChans*surf0m(mif),inst0(mif))
          rad = rad + chanWt0(chan) * surfWt0(mif) * instWt1(mif) * &
            & baseline%values(chan0(chan)+noBslChans*surf0m(mif),inst1(mif))
          rad = rad + chanWt0(chan) * surfWt1(mif) * instWt0(mif) * &
            & baseline%values(chan0(chan)+noBslChans*surf1m(mif),inst0(mif))
          rad = rad + chanWt0(chan) * surfWt1(mif) * instWt1(mif) * &
            & baseline%values(chan0(chan)+noBslChans*surf1m(mif),inst1(mif))
          rad = rad + chanWt1(chan) * surfWt0(mif) * instWt0(mif) * &
            & baseline%values(chan1(chan)+noBslChans*surf0m(mif),inst0(mif))
          rad = rad + chanWt1(chan) * surfWt0(mif) * instWt1(mif) * &
            & baseline%values(chan1(chan)+noBslChans*surf0m(mif),inst1(mif))
          rad = rad + chanWt1(chan) * surfWt1(mif) * instWt0(mif) * &
            & baseline%values(chan1(chan)+noBslChans*surf1m(mif),inst0(mif))
          rad = rad + chanWt1(chan) * surfWt1(mif) * instWt1(mif) * &
            & baseline%values(chan1(chan)+noBslChans*surf1m(mif),inst1(mif))
          radiance%values ( chan + noChans*(mif-1), maf ) = rad
        end do
      end do

      ! --------------------------------------------------------------------
      ! Now compute d[Radiance]/d[Baseline]
      if (present(jacobian) .and. bslInFirst ) then
        instLow = minval(inst0)
        instHi = maxval(inst1)
        allocate ( kBit(radiance%template%instanceLen, &
          & baseline%template%instanceLen, &
          & instLow:instHi), stat=status )
        kBit = 0.0_rp
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate//'kBit' )

        do mif = 1, noMIFs
          mm1 = mif - 1
          do chan = 1, noChans
            row = chan+noChans*mm1
            kBit( row, chan0(chan)+noBslChans*surf0m(mif), inst0(mif) ) = &
              & kBit( row, chan0(chan)+noBslChans*surf0m(mif), inst0(mif) ) + &
              & chanWt0(chan) * surfWt0(mif) * instWt0(mif)
            kBit( row, chan0(chan)+noBslChans*surf0m(mif), inst1(mif) ) = &
              & kBit( row, chan0(chan)+noBslChans*surf0m(mif), inst1(mif) ) + &
              & chanWt0(chan) * surfWt0(mif) * instWt1(mif)
            kBit( row, chan0(chan)+noBslChans*surf1m(mif), inst0(mif) ) = &
              & kBit( row, chan0(chan)+noBslChans*surf1m(mif), inst0(mif) ) + &
              & chanWt0(chan) * surfWt1(mif) * instWt0(mif)
            kBit( row, chan0(chan)+noBslChans*surf1m(mif), inst1(mif) ) = &
              & kBit( row, chan0(chan)+noBslChans*surf1m(mif), inst1(mif) ) + &
              & chanWt0(chan) * surfWt1(mif) * instWt1(mif)
            kBit( row, chan1(chan)+noBslChans*surf0m(mif), inst0(mif) ) = &
              & kBit( row, chan1(chan)+noBslChans*surf0m(mif), inst0(mif) ) + &
              & chanWt1(chan) * surfWt0(mif) * instWt0(mif)
            kBit( row, chan1(chan)+noBslChans*surf0m(mif), inst1(mif) ) = &
              & kBit( row, chan1(chan)+noBslChans*surf0m(mif), inst1(mif) ) + &
              & chanWt1(chan) * surfWt0(mif) * instWt1(mif)
            kBit( row, chan1(chan)+noBslChans*surf1m(mif), inst0(mif) ) = &
              & kBit( row, chan1(chan)+noBslChans*surf1m(mif), inst0(mif) ) + &
              & chanWt1(chan) * surfWt1(mif) * instWt0(mif)
            kBit( row, chan1(chan)+noBslChans*surf1m(mif), inst1(mif) ) = &
              & kBit( row, chan1(chan)+noBslChans*surf1m(mif), inst1(mif) ) + &
              & chanWt1(chan) * surfWt1(mif) * instWt1(mif)
          end do
        end do

        ! Now sparsify and store the blocks
        do instance = lbound(kBit,3), ubound(kBit,3)
          colBlock = FindBlock ( jacobian%col, baseline%index, instance )
          kBit2 => kBit(:,:,instance)
          call Sparsify ( kBit2, jacobian%block(rowBlock,colBlock) )
        end do

        deallocate ( kBit, STAT=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate//'kBit' )
      end if

      ! ---------------------------------------------------------------
      ! Now add some more terms to d[Radiance]/d[ptan] if appropriate
      if (present(jacobian) .and. ptanInFirst) then
        ! Compute the derivative of surfWt[0/1] wrt. ptan
        call Allocate_test ( surfWt0Prime, noMIFs, 'surfWt0Prime', ModuleName )
        call Allocate_test ( surfWt1Prime, noMIFs, 'surfWt1Prime', ModuleName )
      
        where ( surf1 /= surf0 .and. &
          & (ptan%values(:,maf) > baseline%template%surfs(1,1)) .and. &
          & (ptan%values(:,maf) < baseline%template%surfs(baseline%template%noSurfs,1)) )
          surfWt1Prime = 1.0 / &
            & ( baseline%template%surfs(surf1,1) - baseline%template%surfs(surf0,1) ) 
          surfWt0Prime = -1.0 / &
            & ( baseline%template%surfs(surf1,1) - baseline%template%surfs(surf0,1) ) 
        elsewhere
          surfWt1Prime = 0.0
          surfWt0Prime = 0.0
        end where

        if ( any (surfWt0Prime /= 0.0 .or. surfWt1Prime /= 0.0 ) ) then
          colBlock = FindBlock ( jacobian%col, ptan%index, maf )
          jBlock => jacobian%block( rowBlock, colBlock )
          
          ! Create the block as banded if not already done so
          select case (jBlock%kind)
          case (m_absent)
            call CreateBlock ( Jacobian, rowBlock, colBlock, m_banded, &
              & noMIFs*noChans, bandHeight=noChans )
            jBlock%values = 0.0_rp
          case (m_banded)
          case default
            call MLSMessage ( MLSMSG_Error, ModuleName,&
              & 'Wrong matrix type for ptan derivative')
          end select

          ! Now fill in this jacobian block
          do mif = 1, noMIFs
            mm1 = mif - 1
            do chan = 1, noChans
              gradient = 0.0
              gradient = gradient + chanWt0(chan) * surfWt0Prime(mif) * instWt0(mif) * &
                & baseline%values(chan0(chan)+noBslChans*surf0m(mif),inst0(mif))
              gradient = gradient + chanWt0(chan) * surfWt0Prime(mif) * instWt1(mif) * &
                & baseline%values(chan0(chan)+noBslChans*surf0m(mif),inst1(mif))
              gradient = gradient + chanWt0(chan) * surfWt1Prime(mif) * instWt0(mif) * &
                & baseline%values(chan0(chan)+noBslChans*surf1m(mif),inst0(mif))
              gradient = gradient + chanWt0(chan) * surfWt1Prime(mif) * instWt1(mif) * &
                & baseline%values(chan0(chan)+noBslChans*surf1m(mif),inst1(mif))
              gradient = gradient + chanWt1(chan) * surfWt0Prime(mif) * instWt0(mif) * &
                & baseline%values(chan1(chan)+noBslChans*surf0m(mif),inst0(mif))
              gradient = gradient + chanWt1(chan) * surfWt0Prime(mif) * instWt1(mif) * &
                & baseline%values(chan1(chan)+noBslChans*surf0m(mif),inst1(mif))
              gradient = gradient + chanWt1(chan) * surfWt1Prime(mif) * instWt0(mif) * &
                & baseline%values(chan1(chan)+noBslChans*surf1m(mif),inst0(mif))
              gradient = gradient + chanWt1(chan) * surfWt1Prime(mif) * instWt1(mif) * &
                & baseline%values(chan1(chan)+noBslChans*surf1m(mif),inst1(mif))
              jBlock%values ( chan + noChans*mm1, 1 ) = &
                & jBlock%values ( chan + noChans*mm1, 1 ) + gradient
            end do
          end do
        end if                          ! Any nonzero weights
        call Deallocate_test ( surfWt0Prime, 'surfWt0Prime', ModuleName )
        call Deallocate_test ( surfWt1Prime, 'surfWt1Prime', ModuleName )
      end if                             ! Do ptan derivatives

      call Deallocate_test ( inst0, 'inst0', ModuleName )
      call Deallocate_test ( inst1, 'inst1', ModuleName )
      call Deallocate_test ( instWt0, 'instWt0', ModuleName )
      call Deallocate_test ( instWt1, 'instWt1', ModuleName )
      call Deallocate_test ( surf0, 'surf0', ModuleName )
      call Deallocate_test ( surf1, 'surf1', ModuleName )
      call Deallocate_test ( surfWt0, 'surfWt0', ModuleName )
      call Deallocate_test ( surfWt1, 'surfWt1', ModuleName )
      call Deallocate_test ( chan0, 'chan0', ModuleName )
      call Deallocate_test ( chan1, 'chan1', ModuleName )
      call Deallocate_test ( chanWt0, 'chanWt0', ModuleName )
      call Deallocate_test ( chanWt1, 'chanWt1', ModuleName )
    end do                              ! Signal loop

  end subroutine BaselineForwardModel

end module BaselineForwardModel_m
  
! $Log$
! Revision 2.7  2002/01/17 02:16:38  livesey
! Bug, rowBlock wasn't set in some cases
!
! Revision 2.6  2001/10/03 17:46:37  livesey
! Added correction to ptan derivatives
!
! Revision 2.5  2001/10/02 22:22:53  livesey
! Working version
!
! Revision 2.4  2001/10/02 20:37:19  livesey
! Pays attention to do_baseline
!
! Revision 2.3  2001/10/02 20:20:46  livesey
! First linkable version
!
! Revision 2.2  2001/10/02 20:03:58  livesey
! First compilable version
!
! Revision 2.1  2001/10/02 16:52:00  livesey
! Very early version!
!
