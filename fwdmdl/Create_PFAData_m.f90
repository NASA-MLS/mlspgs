! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Create_PFAData_m

  implicit NONE
  private
  public :: Create_PFAData

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  Create_PFAData  -----
  integer function Create_PFAData ( Molecules, Signals, Temperatures, &
    & Pressures, LosVel, WhichLines, Where )

    ! Create PFAData tables for the specified molecules, signals, temperatures
    ! and pressures.  Add them to PFADataBase%PFAData.  Sort PFAData.  Return
    ! the index of the last created PFA datum.

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use DSIMPSON_MODULE, only: SIMPS
    use FilterShapes_m, only: FilterShapes
    use Intrinsic, only: LIT_INDICES, L_NONE
    use L2PC_PFA_STRUCTURES, only: AllocateOneSlabs, DeAllocateOneSlabs, &
      & SLABS_STRUCT
    use MLSCommon, only: RP, R8
    use MLSSignals_m, only: GetNameOfSignal, MatchSignal, MaxSigLen, Signal_T
    use Output_m, only: Output
    use PFADataBase_m, only: AddPFADatumToDatabase, PFAData, PFAData_T, &
      & Sort_PFADataBase
    use Physics, only: h_over_K, SpeedOfLight ! m/s
    use Slabs_SW_m, only: Slabs_Prep_Struct
    use SpectroscopyCatalog_m, only: Catalog, Catalog_t, Line_t, Lines, &
      & MostLines
    use String_Table, only: Display_String
    use Toggles, only: Emit, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use VGridsDatabase, only: VGrid_t

    integer, intent(in) :: Molecules(:)
    type(signal_t), intent(in), target :: Signals(:) ! Derived signals, not from database
    type(vGrid_t), intent(in) :: Temperatures
    type(vGrid_t), intent(in) :: Pressures
    real(rp), intent(in) :: LosVel ! Line-of-sight velocity
    integer, intent(in) :: WhichLines ! 0 => Lines for channel,
                                      ! 1 => Lines for radiometer,
                                      ! 2 => All lines in catalog
    integer, intent(in) :: Where   ! In the parse tree, for error messages
    

    integer :: C, Chan  ! Indices for channels
    integer, pointer :: Channel(:)  ! Index of channel for signal/channel pair
    real(r8) :: DF      ! Spacing in filter bank's frequency grid
    integer :: I        ! Index for signals associated with a line, or lines in catalog
    integer :: L        ! Index for lines in catalog
    logical, pointer :: LINEFLAG(:) ! Use this line
    integer :: M        ! Index for Molecules
    logical, target :: MaxLineFlag(mostLines)
    type(catalog_t) :: MyCatalog
    integer :: N        ! A molecule name
    integer :: NFP      ! Number of points in filter bank's frequency grid
    real(r8) :: Norm    ! Normalization for filter grid
    integer :: NumChannels
    real(rp) :: P       ! Pressure
    type(PFAData_T) :: PFADatum
    logical :: Progress ! Dump signal/molecule
    integer :: PX       ! Index for pressures
    integer :: S        ! Index for signals
    integer :: ShapeInd ! Index of filter shape for signal/sideband/channel
    integer, pointer :: SigInd(:)   ! Index of signal for signal/channel pair
    type(signal_t), pointer :: Signal
    character(maxSigLen) :: SignalName
    logical :: SkipIt   ! No lines or continuum for molecule/signal combination
    type(slabs_struct) :: Slabs
    real(rp) :: T       ! Temperature from temperature grid
    real :: T0, T1, T2  ! for timing
    integer :: TX       ! Index for temperature grid
    type (line_T), pointer :: ThisLine
    real(rp) :: VelCor  ! Velocity correction = 1 - velRel
    real(rp) :: VelRel  ! LosVel/c

    integer, parameter :: NoCat = 1
    integer, parameter :: NoLines = noCat + 1

    if ( toggle(emit) ) & ! set by -f command-line switch
      & call trace_begin ( 'Create_PFAData' )
    progress = index(switches,'pfag') /= 0

    ! Opposite sign convention here from ATBD
    velRel = losVel / speedOfLight ! losVel & speedOfLight both M/s
    velCor = 1.0_rp - velRel

    ! Work out the signal/channel combinations.  Flatten the represenation.
    numChannels = 0
    do s = 1, size(signals)
      numChannels = numChannels + count(signals(s)%channels)
    end do
    nullify ( channel, sigInd )
    call allocate_test ( channel, numChannels, 'Channel', moduleName )
    call allocate_test ( sigInd, numChannels, 'SigInd', moduleName )

    chan = 0
    do s = 1, size(signals)
      do c = 1, size(signals(s)%frequencies)
        if ( signals(s)%channels(c) ) then
          chan = chan + 1
          sigInd(chan) = s
          channel(chan) = c - 1 + lbound(signals(s)%frequencies,1)
        end if
      end do ! c
    end do ! s

    ! Now, for all the channels....
    if ( progress ) then
      call cpu_time ( t0 )
      t1 = t0
    end if
    do c = 1, numChannels
      signal => signals(sigind(c))
      ! Get the filter shape for the signal
      shapeInd = matchSignal ( filterShapes%signal, signal, &
        & sideband=signal%sideband, channel=channel(c) )
      pfaDatum%filterFile = filterShapes(shapeInd)%file
      nfp = size(filterShapes(shapeInd)%filterGrid)
      df = filterShapes(shapeInd)%filterGrid(2) - filterShapes(shapeInd)%filterGrid(1)
      ! Compute integral of filter shape, for normalization.  Should be 1.0,
      ! but maybe the input file didn't get normalized....
      call simps ( filterShapes(shapeInd)%filterShape, df, nfp, norm )

      ! Now, for all the molecules....
      do m = 1, size(molecules)
        n = molecules(m)

        ! Create an empty PFA Datum
        pfaDatum%name = 0
        pfaDatum%molecule = n
        call getNameOfSignal ( signal, pfaDatum%signal, channel=channel(c) )
        pfaDatum%signalIndex = signal%index
        pfaDatum%theSignal = signal
        nullify ( pfaDatum%theSignal%channels )
        call allocate_test ( pfaDatum%theSignal%channels, channel(c), &
          & 'PFADatum%theSignal%channels', moduleName, lowBound=channel(c) )
        pfaDatum%theSignal%channels = .true.
        pfaDatum%tGrid = temperatures
        pfaDatum%vGrid = pressures
        pfaDatum%vel_Rel = velRel

        px = pfaDatum%vGrid%noSurfs
        tx = pfaDatum%tGrid%noSurfs
        nullify ( pfaDatum%absorption, pfaDatum%dAbsDnc, pfaDatum%dAbsDnu, pfaDatum%dAbsDwc )
        call allocate_test ( pfaDatum%absorption, tx, px, 'pfaDatum%absorption', moduleName )
        call allocate_test ( pfaDatum%dAbsDnc,    tx, px, 'pfaDatum%dAbsDnc',    moduleName )
        call allocate_test ( pfaDatum%dAbsDnu,    tx, px, 'pfaDatum%dAbsDnu',    moduleName )
        call allocate_test ( pfaDatum%dAbsDwc,    tx, px, 'pfaDatum%dAbsDwc',    moduleName )

        ! Work out the required spectroscopy
        myCatalog = catalog(n)
        if ( myCatalog%molecule == l_none ) then
          call announce_error ( where, noCat, more=lit_indices(n) )
          cycle
        end if
        call work_out_spectroscopy
        if ( skipIt ) then
          if ( progress ) call cpu_time ( t1 )
          cycle
        end if

        ! Then for all the pressures and temperatures....
        do tx = 1, temperatures%noSurfs
          t = exp(temperatures%surfs(tx,1))
          do px = 1, pressures%noSurfs
            p = 10.0_rp**(-pressures%surfs(px,1))
            call allocateOneSlabs ( slabs, myCatalog, .false. )
            call slabs_prep_struct ( t, p, myCatalog, velCor, .false., slabs )
            call get_beta_etc
            call deallocateOneSlabs ( slabs, moduleName )
          end do ! px
        end do ! tx

        ! Put it away
        create_PFAData = AddPFADatumToDatabase ( pfaData, pfaDatum )

        if ( progress ) then
          call output ( 'Created PFA for ' )
          call display_string ( lit_indices(n) )
          call output ( ' / ' )
          call output ( trim(pfaDatum%signal) )
          call cpu_time ( t2 )
          call output ( t2-t1, before=' using ', after=' seconds', &
            & format='(f0.2)', advance='yes' )
          t1 = t2
        end if
      end do ! m
    end do ! c

    call deallocate_test ( channel, 'Channel', moduleName )
    call deallocate_test ( sigInd, 'SigInd', moduleName )

    call sort_PFADataBase

    if ( progress ) then
      call cpu_time ( t2 )
      call output ( t2-t0, before='Total CPU time = ', advance='yes' )
    end if

    if ( toggle(emit) ) & ! set by -f command-line switch
      & call trace_end ( 'Create_PFAData' )

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String, More )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      use String_Table, only: Display_String
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String
      integer, intent(in), optional :: More
      call startErrorMessage ( where )
      select case ( what )
      case ( noCat )
        call output ( 'No catalog for ' )
        call display_string ( more, advance='yes' )
      case ( noLines )
        call output ( 'No lines or continuum for ' )
        call display_string ( more )
        call output ( ' / ' )
        call output ( trim(string), advance='yes' )
      end select
    end subroutine Announce_Error

    ! .............................................  Get_Beta_Etc  .....
    subroutine Get_Beta_Etc
      ! Get Beta, dBeta_dw, dBeta_dn and dBeta_dv.  Don't need dBeta_dT.
      ! Frequency average them.
      use Get_Beta_Path_m, only: Create_Beta
      real(r8) :: Avg, dAvg
      real(rp), dimension(size(filterShapes(shapeInd)%filterGrid)) :: &
        & Beta, dBeta_dw, dBeta_dn, dBeta_dv
      integer :: F    ! Index for frequencies
      real(r8) :: FRQ ! Frequency from filter grid
      real(rp), parameter :: h_over_2K = 0.5 * h_over_K
      real(r8) :: Temp(size(filterShapes(shapeInd)%filterGrid))
      ! Compute Beta and its derivatives
      do f = 1, nfp
        frq = filterShapes(shapeInd)%filterGrid(f)
        beta(f) = 0.0
        call create_beta ( p, T, frq, slabs, &
          & real(tanh(h_over_2K * frq / T),rp), beta(f), noPolarized=.false., &
          & dBeta_dw=dBeta_dw(f), dBeta_dn=dBeta_dn(f), dBeta_dv=dBeta_dv(f) )
      end do ! f
      ! Average.  Assumes filter grid's frequencies are evenly spaced.
      temp = beta * filterShapes(shapeInd)%filterShape
      call simps ( temp, df, nfp, avg )
      pfaDatum%Absorption(tx,px) = log( avg / norm )
      temp = dBeta_dw * filterShapes(shapeInd)%filterShape
      call simps ( temp, df, nfp, dAvg )  ! normalization cancels for derivs
      pfaDatum%dAbsDwc(tx,px) = dAvg / avg ! d ln beta / d w = 1 / beta d beta / d w
      temp = dBeta_dn * filterShapes(shapeInd)%filterShape
      call simps ( temp, df, nfp, dAvg )
      pfaDatum%dAbsDnc(tx,px) = dAvg / avg ! d ln beta / d n = 1 / beta d beta / d n
      temp = dBeta_dv * filterShapes(shapeInd)%filterShape
      call simps ( temp, df, nfp, dAvg )
      pfaDatum%dAbsDnu(tx,px) = dAvg / avg ! d ln beta / d v = 1 / beta d beta / d v
    end subroutine Get_Beta_Etc

    ! ....................................  Work_Out_Spectroscopy  .....
    subroutine Work_Out_Spectroscopy
      use MLSSignals_m, only: GetRadiometerFromSignal
      skipIt = .false. ! Assume there will be lines and/or continuum
      ! Don't deallocate lines by mistake -- myCatalog is a shallow copy
      nullify ( myCatalog%lines )
      if ( associated ( catalog(n)%lines ) ) then
        ! Subset the lines according to the signal
        lineFlag => MaxLineFlag(:size(catalog(n)%lines))
        if ( whichLines > 1 ) then
          lineFlag = .true.
        else
          do l = 1, size ( catalog(n)%lines )
            lineFlag(l) = .false.
            thisLine => lines(catalog(n)%lines(l))
            if ( associated(thisLine%signals) ) then
              do i = 1, size(thisLine%signals)
                if ( whichLines > 0 .and. &
                    & getRadiometerFromSignal(thisLine%signals(i)) == &
                    & signal%radiometer .or. &
                  & thisLine%signals(i) == signal%index .and. &
                  &  ( thisLine%sidebands(i) * signal%sideband == 0 .or. &
                  &    thisline%sidebands(i) == signal%sideband ) ) then
                  lineFlag(l) = .true.
                  exit
                end if
              end do ! i
            end if ! associated(thisLine%signals)
          end do ! l
        end if
        ! Check we have at least one line for specie.  Allocate lines if so.
        l = count(lineFlag)
        if ( l == 0 .and. all(myCatalog%continuum == 0) ) then
          call getNameOfSignal ( signal, signalName )
          call announce_error ( where, noLines, signalName, lit_indices(n) )
          skipIt = .true.
          return
        end if
        call allocate_test ( myCatalog%lines, l, 'myCatalog%lines', moduleName )
        myCatalog%lines = pack ( catalog(n)%lines, lineFlag )

      else

        ! No lines for this specie.  However, its continuum is still valid 
        ! so don't set it to empty.
        ! Don't bother checking that continuum /= 0 as otherwise then
        ! presumably having no continuum and no lines it wouldn't be in the
        ! catalog!
        call allocate_test ( myCatalog%lines, 0, 'myCatalog%lines(0)', moduleName )
      end if
    end subroutine Work_Out_Spectroscopy
  end function Create_PFAData

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Create_PFAData_m

! $Log$
! Revision 2.5  2005/03/16 23:59:56  vsnyder
! Add allLinesForRadiometer and allLinesInCatalog to makePFA
!
! Revision 2.4  2005/01/27 21:20:08  vsnyder
! Remove nonscalar molecule
!
! Revision 2.3  2005/01/12 03:17:26  vsnyder
! Set the channel correctly
!
! Revision 2.2  2004/12/31 02:41:01  vsnyder
! Create PFA data
!
! Revision 2.1  2004/12/13 23:54:46  vsnyder
! Initial commit
!
