! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Create_PFAData_m

  implicit NONE
  private
  public :: Create_PFAData

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  Create_PFAData  -----
  integer function Create_PFAData ( Molecules, Signals, Temperatures, &
    & Pressures, LosVel, WhichLines, Oversample, Where )

    ! Create PFAData tables for the specified molecules, signals, temperatures
    ! and pressures.  Tables for DACS are calculated using unpolarized betas
    ! (we have no provision for magnetic field in the PFA tables).  Add the
    ! tables to PFADataBase%PFAData.  Return the index of the last created PFA
    ! datum.

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSNUMERICS, only: SIMPS => SIMPSONSSUB
    use DUMP_0, only: DUMP
    use FILTERSHAPES_M, only: DACSFILTERSHAPES, FILTERSHAPES, FILTERSHAPE_T
    use INTRINSIC, only: LIT_INDICES, L_NONE
    use MLSKINDS, only: RP, R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSSIGNALS_M, only: GETNAMEOFSIGNAL, MATCHSIGNAL, MAXSIGLEN, SIGNAL_T
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use MORETREE, only: GETSTRINGINDEXFROMSTRING
    use OUTPUT_M, only: OUTPUT
    use PFADATABASE_M, only: ADDPFADATUMTODATABASE, HOOKTABLETOFINDPFA, &
      & PFADATA, PFADATA_T
    use PHYSICS, only: H_OVER_K, SPEEDOFLIGHT ! M/S
    use SLABS_SW_M, only: ALLOCATEONESLABS, DEALLOCATEONESLABS, &
      & SLABS_STRUCT, SLABS_PREP_STRUCT
    use SPECTROSCOPYCATALOG_M, only: CATALOG, CATALOG_T, LINE_T, LINES, &
      & MOSTLINES, SPECTROSCOPYFILE
    use STRING_TABLE, only: DISPLAY_STRING
    use TOGGLES, only: EMIT, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VGRIDSDATABASE, only: VGRID_T

    type(filterShape_t), pointer :: Filters(:)
    integer, intent(in) :: Molecules(:)
    type(filterShape_t), pointer :: MyFilter
    type(signal_t), intent(in), target :: Signals(:) ! Derived signals, not from database
    type(vGrid_t), intent(in) :: Temperatures
    type(vGrid_t), intent(in) :: Pressures
    real(rp), intent(in) :: LosVel ! Line-of-sight velocity
    integer, intent(in) :: WhichLines ! 0 => Lines for channel,
                                      ! 1 => Lines for radiometer,
                                      ! 2 => All lines in catalog
    integer, intent(in) :: Oversample ! How much to oversample the filter grid
    integer, intent(in) :: Where   ! In the parse tree, for error messages
    

    integer :: C, Chan  ! Indices for channels
    integer, pointer :: Channel(:)  ! Index of channel for signal/channel pair
    real(r8) :: DF      ! Spacing in filter bank's frequency grid
    integer :: DumpIt   ! Dump Beta, dBetaD... if nonzero.  Stop after first
                        !   one if > 1.
    logical :: Error
    integer :: I        ! Index for signals associated with a line, or lines in catalog
    integer :: Ix       ! Position in PFAData of a sought PFADatum if nonzero
    logical, pointer :: LINEFLAG(:) ! Use this line
    integer :: M        ! Index for Molecules
    logical, target :: MaxLineFlag(mostLines)
    integer :: Me = -1  ! String index for trace
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
    character(len=maxSigLen) :: SignalText
    logical :: SkipIt   ! No lines or continuum for molecule/signal combination
    type(slabs_struct) :: Slabs
    real(rp) :: T       ! Temperature from temperature grid
    real :: T0, T1, T2  ! for timing
    integer :: TX       ! Index for temperature grid
    real(rp) :: VelCor  ! Velocity correction = 1 - velRel
    real(rp) :: VelRel  ! LosVel/c

    integer, parameter :: DupSignals = 1
    integer, parameter :: NoCat = dupSignals + 1
    integer, parameter :: NoFilter = noCat + 1
    integer, parameter :: NoLines = noFilter + 1


    call trace_begin ( me, 'Create_PFAData', &
      & cond=toggle(emit) ) ! set by -f command-line switch
    progress = switchDetail(switches,'pfag') > -1
    dumpIt =  switchDetail(switches,'pfab')

    error = .false.

    ! Opposite sign convention here from ATBD
    velRel = losVel / speedOfLight ! losVel & speedOfLight both M/s
    velCor = 1.0_rp - velRel

    ! Work out the signal/channel combinations.  Flatten the representation.
    ! Check for duplicate signals
    numChannels = 0
    do s = 1, size(signals)
      numChannels = numChannels + count(signals(s)%channels)
      do i = 1, s-1
        if ( matchSignal(signals(s),signals(i)) > 0 ) then
          call getNameOfSignal ( signals(s), signalText )
          call announce_error ( where, dupSignals, signalText )
          error = .true.
        end if
      end do
    end do
    nullify ( channel, sigInd )
    if ( error ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, 'Duplicate signals' )

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
    pfaDatum%name = 0
    do c = 1, numChannels
      signal => signals(sigind(c))

      ! Create an empty PFA Datum
      call getNameOfSignal ( signal, signalText, channel=channel(c) )
      pfaDatum%channel = channel(c)
      pfaDatum%signal = getStringIndexFromString ( trim(signalText) )
      pfaDatum%signalIndex = signal%index
      pfaDatum%spectroscopyFile = spectroscopyFile
      pfaDatum%theSignal = signal
      pfaDatum%tGrid = temperatures
      pfaDatum%vel_Rel = velRel
      pfaDatum%vGrid = pressures
      pfaDatum%whichLines = whichLines
      ! Get the filter shape for the signal
      nullify ( filters )
      if ( signal%dacs ) then
        if (  associated(DACSFilterShapes)) filters => DACSFilterShapes%filter
      else
        filters => filterShapes
      end if
      if ( .not. associated(filters) ) &
        & call announce_error ( where, noFilter, signalText )
      shapeInd = matchSignal ( filters%signal, signal, &
        & sideband=signal%sideband, channel=channel(c) )
      if ( shapeInd == 0 ) then
        call announce_error ( where, noFilter, signalText )
        cycle
      end if
      myFilter => filters(shapeInd)
      pfaDatum%filterFile = myFilter%file
      nfp = size(myFilter%filterGrid)
      df = myFilter%filterGrid(2) - myFilter%filterGrid(1)
      ! Compute integral of filter shape, for normalization.  Should be 1.0,
      ! but maybe the input file didn't get normalized....
      call simps ( myFilter%filterShape, df, nfp, norm )
      df = df / oversample

      ! Now, for all the molecules....
      do m = 1, size(molecules)
        n = molecules(m)
        pfaDatum%molecule = n
        ! The channels field is allocated here so each one in the
        ! database will have a separate one, even if they are the same.
        ! Otherwise, when it comes time to destroy them, the first one
        ! will work, and the next one will fail with a dangling pointer.
        nullify ( pfaDatum%theSignal%channels )
        call allocate_test ( pfaDatum%theSignal%channels, channel(c), &
          & 'PFADatum%theSignal%channels', moduleName, lowBound=channel(c) )
        pfaDatum%theSignal%channels = .true.
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
          ! temperatures%surfs is really log temperature
          t = exp(temperatures%surfs(tx,1))
          do px = 1, pressures%noSurfs
            p = 10.0_rp**(-pressures%surfs(px,1))
            call allocateOneSlabs ( slabs, myCatalog, moduleName, .false. )
            call slabs_prep_struct ( t, p, myCatalog, velCor, .false., slabs )
            call get_beta_etc
            call deallocateOneSlabs ( slabs, moduleName )
          end do ! px
        end do ! tx

        ! Look for it
        ix = HookTableToFindPFA ( 0, 0, PFADatum, 0, replace=.true. )

        ! Put it away
        if ( ix /= 0 ) then
          PFAData(ix) = PFADatum
        else
          create_PFAData = AddPFADatumToDatabase ( PFAData, PFADatum )
          call deallocate_test ( myCatalog%lines, 'myCatalog%lines', moduleName )
          ix = HookTableToFindPFA ( 0, 0, PFADatum, create_PFAData )
        end if

        if ( progress .or. dumpIt >= 0 ) then
          call output ( 'Created PFA for ' )
          call display_string ( lit_indices(n) )
          call output ( ' / ' )
          call output ( trim(signalText) )
          call cpu_time ( t2 )
          call output ( t2-t1, before=' using ', after=' seconds', &
            & format='(f0.2)', advance='yes' )
          t1 = t2
        end if

        if ( dumpIt >= 0 ) then
          call dump ( pfaDatum%Absorption, name='Absorption' )
          call dump ( pfaDatum%dAbsDwc, name='DAbsDwc' )
          call dump ( pfaDatum%dAbsDnc, name='DAbsDnc' )
          call dump ( pfaDatum%dAbsDnu, name='DAbsDnu' )
          if ( dumpIt > 1 ) stop
        end if
      end do ! m
    end do ! c
    if ( error ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Errors while hooking up and finding PFA tables' )

    call deallocate_test ( channel, 'Channel', moduleName )
    call deallocate_test ( sigInd, 'SigInd', moduleName )

    if ( progress ) then
      call cpu_time ( t2 )
      call output ( t2-t0, before='Total CPU time for CreatePFA = ', advance='yes' )
    end if

    call trace_end ( 'Create_PFAData', cond=toggle(emit) )

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String, More )
      use MORETREE, only: STARTERRORMESSAGE
      use OUTPUT_M, only: OUTPUT
      use STRING_TABLE, only: DISPLAY_STRING
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String
      integer, intent(in), optional :: More
      call startErrorMessage ( where )
      select case ( what )
      case ( dupSignals )
        call output ( 'Signal ' )
        call output ( trim(string) )
        call output (' duplicates a previously-specified signal.', advance='yes' )
      case ( noCat )
        call output ( 'No catalog for ' )
        call display_string ( more, advance='yes' )
      case ( noFilter )
        call output ( 'No filter for ' )
        call output ( trim(string), advance='yes' )
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
      real(rp), dimension((nfp-1)*oversample+1) :: &
        & Beta, dBeta_dw, dBeta_dn, dBeta_dv, Frqs, Shapes
      integer :: F    ! Index for frequencies
      real(r8) :: FRQ ! Frequency from filter grid
      real(rp), parameter :: h_over_2K = 0.5 * h_over_K
      integer :: I    ! Index for filling oversamples
      integer :: NQP  ! Number of quadrature points = (nfp-1)*oversample+1
      real(r8) :: Temp(size(Shapes))
      ! Oversample the filter grid using linear interpolation
      nqp = size(frqs)
      frqs(1) = myFilter%filterGrid(1)
      shapes(1) = myFilter%filterShape(1)
      do f = 1+oversample, nqp, oversample
        frqs(f) = myFilter%filterGrid((f+oversample-1)/oversample)
        shapes(f) = myFilter%filterShape((f+oversample-1)/oversample)
        do i = 1, oversample-1
          frqs(f-oversample+i) = frqs(f-oversample) + &
            & real(i)/oversample * (frqs(f)-frqs(f-oversample))
          shapes(f-oversample+i) = shapes(f-oversample) + &
            & real(i)/oversample * (shapes(f)-shapes(f-oversample))
        end do ! i
      end do ! f
      ! Compute Beta and its derivatives
      do f = 1, nqp
        frq = frqs(f)
        beta(f) = 0.0
        call create_beta ( p, T, frq, signals(1)%lo, slabs, &
          & real(tanh(h_over_2K * frq / T),rp), beta(f), noPolarized=.false., &
          & dBeta_dw=dBeta_dw(f), dBeta_dn=dBeta_dn(f), dBeta_dv=dBeta_dv(f) )
      end do ! f
      ! Average.  Assumes filter grid's frequencies are evenly spaced.
      temp = beta * shapes
      call simps ( temp, df, nqp, avg )
      pfaDatum%Absorption(tx,px) = log( avg / norm )
      temp = dBeta_dw * shapes
      call simps ( temp, df, nqp, dAvg )  ! normalization cancels for derivs
      pfaDatum%dAbsDwc(tx,px) = dAvg / avg ! d ln beta / d w = 1 / beta d beta / d w
      temp = dBeta_dn * shapes
      call simps ( temp, df, nqp, dAvg )
      pfaDatum%dAbsDnc(tx,px) = dAvg / avg ! d ln beta / d n = 1 / beta d beta / d n
      temp = dBeta_dv * shapes
      call simps ( temp, df, nqp, dAvg )
      pfaDatum%dAbsDnu(tx,px) = dAvg / avg ! d ln beta / d v = 1 / beta d beta / d v
    end subroutine Get_Beta_Etc

    ! ....................................  Work_Out_Spectroscopy  .....
    subroutine Work_Out_Spectroscopy
      use MLSSignals_m, only: GETRADIOMETERFROMSIGNAL
      integer :: I, L
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
            associate ( thisLine => lines(catalog(n)%lines(l)) )
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
            end associate
          end do ! l
        end if
        ! Check we have at least one line for specie.  Allocate lines if so.
        l = count(lineFlag)
        if ( l == 0 .and. all(myCatalog%continuum == 0) ) then
          call announce_error ( where, noLines, signalText, lit_indices(n) )
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Create_PFAData_m

! $Log$
! Revision 2.30  2018/08/06 19:56:48  vsnyder
! Use ASSOCIATE construct to avoid necessity for Lines database to have the
! TARGET attribute.  Some cannonball polishing.
!
! Revision 2.29  2017/12/07 02:42:49  vsnyder
! Don't use host-associated DO indices; make them local
!
! Revision 2.28  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.27  2013/07/26 22:19:04  vsnyder
! Fiddle with dump switches
!
! Revision 2.26  2011/08/26 00:31:39  pwagner
! CSpline and Hunt now USE MLSNumerics
!
! Revision 2.25  2011/05/09 17:44:26  pwagner
! Converted to using switchDetail
!
! Revision 2.24  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.23  2008/10/03 16:26:32  livesey
! Pushed lo down to support EXTINCTIONV2
!
! Revision 2.22  2008/09/04 19:59:32  vsnyder
! Add PRINT statement in not_used_here
!
! Revision 2.21  2006/07/29 03:01:41  vsnyder
! Send module name to AllocateOneSlabs
!
! Revision 2.20  2006/07/10 22:26:20  vsnyder
! Integrate entire grid when oversampling
!
! Revision 2.19  2006/07/08 01:14:17  vsnyder
! Don't create junk at end of oversampled arrays
!
! Revision 2.18  2006/06/15 20:39:29  vsnyder
! Add oversampling
!
! Revision 2.17  2006/04/26 00:39:26  vsnyder
! Need either ordinary or DACS filters
!
! Revision 2.16  2006/04/25 23:25:36  vsnyder
! Revise DACS filter shape data structure
!
! Revision 2.15  2006/04/21 22:24:19  vsnyder
! Stuff for updating PFA
!
! Revision 2.14  2006/01/26 03:08:05  vsnyder
! Accumulate all errors before crashing, detect duplicate signals early
!
! Revision 2.13  2005/06/09 02:34:15  vsnyder
! Move stuff from l2pc_pfa_structures to slabs_sw_m
!
! Revision 2.12  2005/06/03 01:58:53  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! Revise PFA data structures.
!
! Revision 2.11  2005/05/24 01:54:29  vsnyder
! Delete unused symbols
!
! Revision 2.10  2005/05/12 20:49:28  livesey
! Bug fix (hopefully), fill in channel information in PFADatum on
! creation.
!
! Revision 2.9  2005/05/02 23:02:18  vsnyder
! Use string index of signal's text
!
! Revision 2.8  2005/04/04 19:52:46  vsnyder
! Hoist some loop-invariant stuff
!
! Revision 2.7  2005/03/28 20:22:59  vsnyder
! Add more progress dumps
!
! Revision 2.6  2005/03/17 01:32:26  vsnyder
! Put spectroscopy file's string index in PFAData structure
!
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
