! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module AntennaPatterns_m

  ! Read the antenna patterns file.

  use MLSKINDS, only: R8
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
  use MLSSIGNALS_M, only: MAXSIGLEN, SIGNALS, SIGNAL_T, DISPLAYSIGNALNAME

  implicit NONE

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Antenna_Patterns_File, Read_Antenna_Patterns_File
  public :: Close_Antenna_Patterns_File, Destroy_Ant_Patterns_Database
  public :: Dump_Antenna_Patterns_Database

  type, public :: AntennaPattern_T
    real(r8) :: Lambda
!     ! FFTW_*_Plan are actually places for FFTW to stash C pointers
!     real(r8) :: FFTW_Forward_Plan, FFTW_Backward_Plan
!     ! Aaap, DD1aap and D2aap are produced by Math77 DRFT1.
!     ! The first two elements are the real parts of the first and last
!     ! Fourier coefficients, for which the imaginary parts are zero.
!     ! Remaining pairs of elements are real and imaginary parts of
!     ! subsequent coefficients.  So (3:4) are the real and imaginary
!     ! parts of the second coefficient, (5:6) for the third, etc.
    real(r8), dimension(:), pointer :: Aaap => NULL()
    real(r8), dimension(:), pointer :: D1aap => NULL()
    real(r8), dimension(:), pointer :: D2aap => NULL()
    type (Signal_T), pointer, dimension(:) :: Signals => NULL()
  end type AntennaPattern_T

  ! The antenna pattern database:
  type(AntennaPattern_T), dimension(:), pointer, save, public :: &
    & AntennaPatterns => NULL()

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Antenna_Patterns_File  -----
  subroutine Open_Antenna_Patterns_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the antenna pattern file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    integer :: Status

    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open antenna pattern file " // Filename )
  end subroutine Open_Antenna_Patterns_File

  ! ------------------------------------  Read_Antenna_Patterns_File  -----
  subroutine Read_Antenna_Patterns_File ( Lun, Where )
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST, &
      & Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use MACHINE, only: IO_ERROR
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use PARSE_SIGNAL_M, only: PARSE_SIGNAL
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use CONSTANTS, only: PI

    integer, intent(in) :: Lun               ! Logical unit number to read it
    integer, intent(in) :: Where             ! In the L2CF tree, for tracing

!     ! Parameters for FFTW
!     integer, parameter :: FFTW_Forward = -1, FFTW_Backward = 1
!     integer, parameter :: FFTW_Measure = 1
! 
!     ! Interface for RFFTW_f77_create_plan
!     interface
!       subroutine RFFTW_f77_create_plan ( Plan, NT, Dir, Flags )
!         use MLSCommon, only: R8
!         real(r8), intent(out) :: Plan ! Actually a C pointer
!         integer , intent(in) :: NT    ! Number of points
!         integer , intent(in) :: Dir   ! Direction, forward or backward
!         integer , intent(in) :: Flags ! FFTW input flags
!       end subroutine RFFTW_f77_create_plan
!     end interface

    real(r8), parameter :: Pi2 = 2.0_r8 * Pi

    integer(c_intptr_t) :: Addr         ! For tracing
    logical, dimension(:), pointer :: CHANNELS ! From Parse Signal
    integer :: DataBaseSize                  ! How many antenna patterns?
    integer :: HowManyPoints(3*size(signals))  ! for each pattern
    integer :: HowManySignals(3*size(signals)) ! for each pattern
    integer :: HowManySignalLines(3*size(signals)) ! for each pattern
    integer :: I, J, K, L                    ! Loop inductors, subscripts
    integer :: Me = -1                       ! String index for trace
    integer :: power2                        ! Power 2 size of antenna pattern
    real(r8) :: Lambda
    real(r8) :: LambdaX2Pi                   ! 2 * Pi * Lambda
    real(r8) :: Q                            ! Factor used to scale derivatives
    character(len=MaxSigLen) :: SigName      ! Signal Name
    integer :: SIDEBAND                      ! From parse_signal
    integer :: Status                        ! From read or allocate
    integer :: DummyInt                      ! Dummy integer to read
    integer, dimension(:), pointer :: SIGINDS ! From parseSignal
    integer :: SignalCount
    integer, pointer, dimension(:) :: Signal_Indices => NULL() ! From Parse_Signal, q.v.
    !                                          It's never allocated because of
    !                                          the "only_count_em" argument
    real(r8) :: V(2)                         ! To read a line from the file

    call trace_begin ( me, "Read_Antenna_Patterns_File", where, cond=toggle(gen) )

    if ( associated(AntennaPatterns) ) call destroy_ant_patterns_database


    ! First, read through the file and count how much stuff is there.
    read ( lun, '(a)', end=98, err=99, iostat=status ) sigName
    dataBaseSize = 0
outer1: do
      dataBaseSize = dataBaseSize + 1
      if ( dataBaseSize > size(howManySignals) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & "More patterns in the file than I have space for" )
      howManySignals(dataBaseSize) = 0
      howManySignalLines(dataBaseSize) = 0
      do ! Count how many signals there are
        sigName = adjustl(sigName)
        if ( verify(sigName(1:1), '0123456789.+-') == 0 ) exit ! a number
        call parse_signal ( sigName, signal_indices, &
          & onlyCountEm=signalCount )
        if ( signalCount == 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & trim(sigName) // " is not a valid signal." )
        howManySignalLines(dataBaseSize) = howManySignalLines(dataBaseSize) + 1
        howManySignals(dataBaseSize) = howManySignals(dataBaseSize) + signalCount
        read ( lun, '(a)', end=98, err=99, iostat=status ) sigName
      end do
      read ( sigName, *, err=99, iostat=status ) dummyInt, lambda
      howManyPoints(dataBaseSize) = 0
      do ! Count how many points in the pattern
        read ( lun, '(a)', err=99, iostat=status ) sigName
        if ( status < 0 ) exit outer1
        sigName = adjustl(sigName)
        if ( verify(sigName(1:1), '0123456789.+-') /= 0 ) exit ! not a number
        read ( sigName, *, err=99, iostat=status ) v ! Check the numbers
        howManyPoints(dataBaseSize) = howManyPoints(dataBaseSize) + 1
      end do
    end do outer1

    ! Now we know how big the database is, and each part of it.
    allocate ( antennaPatterns(dataBaseSize), stat=status )
    addr = 0
    if ( status == 0 .and. dataBaseSize > 0 ) &
      & addr = transfer(c_loc(antennaPatterns(1)), addr)
    call test_allocate ( status, moduleName, "AntennaPatterns", &
      & uBounds = dataBaseSize, elementSize = storage_size(antennaPatterns) / 8, &
      & address=addr )
    rewind ( lun )
    do i = 1, dataBaseSize
      allocate  ( antennaPatterns(i)%signals(howManySignals(i)), stat=status )
      addr = 0
      if ( status == 0 .and. howManySignals(i) > 0 ) &
        & addr = transfer(c_loc(antennaPatterns(i)%signals(1)), addr)
      call test_allocate ( status, moduleName, "AntennaPatterns%signals", &
        & uBounds = howManySignals(i), &
        & elementSize = storage_size(antennaPatterns(i)%signals) / 8, address=addr )

      power2 = 2**ceiling( log10( real(howManyPoints(i)) ) / log10(2.0) )
! billsdebug
      power2 = 2**11

      call Allocate_Test ( antennaPatterns(i)%aaap, 2*power2, &
        & "AntennaPatterns(?)%Aaap", moduleName )
      call Allocate_Test ( antennaPatterns(i)%d1aap, 2*power2, &
        & "AntennaPatterns(?)%D1aap", moduleName )
      call Allocate_Test ( antennaPatterns(i)%d2aap, 2*power2, &
        & "AntennaPatterns(?)%D2aap", moduleName )

      ! Fill with zeros, then we'll read the first part
      antennaPatterns(i)%aaap = 0.0
      antennaPatterns(i)%d1aap = 0.0
      antennaPatterns(i)%d2aap = 0.0

!       ! Find FFTW plans for the current size, or create new ones
!       do j = 1, i-1
!         if ( size(antennaPatterns(i)%aaap) == size(antennaPatterns(j)%aaap) ) exit
!       end do
! 
!       if ( j < i ) then
!         antennaPatterns(i)%FFTW_Forward_Plan = antennaPatterns(j)%FFTW_Forward_Plan
!         antennaPatterns(i)%FFTW_Backward_Plan = antennaPatterns(j)%FFTW_Backward_Plan
!       else
!         call rfftw_f77_create_plan ( antennaPatterns(i)%FFTW_Forward_Plan, &
!           & size(antennaPatterns(i)%aaap), FFTW_Forward, FFTW_Measure )
!         call rfftw_f77_create_plan ( antennaPatterns(i)%FFTW_Backward_Plan, &
!           & size(antennaPatterns(i)%aaap), FFTW_Backward, FFTW_Measure )
!       end if

      k = 1
      do j = 1, howManySignalLines(i)
        read ( lun, '(a)', err=99, iostat=status ) SigName
        nullify ( channels, sigInds )
        call Parse_Signal ( SigName, sigInds, sideband=sideband, channels=channels )
        antennaPatterns(i)%signals(k:k+size(sigInds)-1) = &
          & signals(sigInds)
        antennaPatterns(i)%signals(k:k+size(sigInds)-1)%sideband = sideband
        if (associated(channels)) then
          do l = 1, size(sigInds)
            ! Now nullify the channels we have.  Don't hose the main database!
            nullify ( antennaPatterns(i)%signals(k+l-1)%channels )
            call allocate_test ( antennaPatterns(i)%signals(k+l-1)%channels, &
              & size(antennaPatterns(i)%signals(k+l-1)%frequencies), &
              & 'antennaPatterns(?)%signals(?)%channels', &
              & moduleName )
            antennaPatterns(i)%signals(k+l-1)%channels = .false.
            antennaPatterns(i)%signals(k+l-1)%channels ( &
              & lbound(channels,1):ubound(channels,1) ) = channels
          end do
        end if
        k = k + size(sigInds)
        call deallocate_test ( sigInds, 'sigInds', ModuleName )
        call deallocate_test ( channels, 'channels', ModuleName )
      end do ! j
      read ( lun, *, err=99, iostat=status ) dummyInt, lambda
      antennaPatterns(i)%lambda = lambda
      lambdaX2Pi = pi2 * lambda
      do j = 1, howManyPoints(i)
        read ( lun, *, err=99, iostat=status ) v
        k = 2 * j - 1
        antennaPatterns(i)%aaap(k:k+1) = v(1:2)

        ! First derivative field:     i*Q * F(S), i = Sqrt(-1)

        q = (j-1) * lambdaX2Pi
        antennaPatterns(i)%d1aap(k)    = -v(2) * q
        antennaPatterns(i)%d1aap(k+1)  =  v(1) * q

        ! Second derivative field:    (i*Q)**2 * F(S), i = Sqrt(-1)

        q = -q * q
        antennaPatterns(i)%d2aap(k:k+1)  =  v(1:2) * q
      end do ! j
    end do ! i

    if ( switchDetail(switches,'ant') > -1 ) call dump_Antenna_Patterns_database
    call trace_end ( "Read_Antenna_Patterns_File", cond=toggle(gen) )

    return
 98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
 99 call io_error ( "While reading the antenna pattern file", status )
    call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )
  end subroutine Read_Antenna_Patterns_File

  ! -----------------------------------  Close_Antenna_Patterns_File  -----
  subroutine Close_Antenna_Patterns_File ( Lun )
    integer, intent(in) :: lun
    close ( lun )
  end subroutine Close_Antenna_Patterns_File

  ! ---------------------------------  Destroy_Ant_Patterns_Database  -----
  subroutine Destroy_Ant_Patterns_Database
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

!     interface
!       subroutine RFFTW_f77_destroy_plan ( Plan )
!         use MLSCommon, only: R8
!         real(r8), intent(inout) :: Plan
!       end subroutine RFFTW_f77_destroy_plan
!     end interface

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, J, S, Status

    if (.not. associated(AntennaPatterns) ) return

!     ! Destroy FFTW plans
!     do i = 1, size(AntennaPatterns)
!       do j = 1, i-1
!         if ( size(antennaPatterns(i)%aaap) == size(antennaPatterns(j)%aaap) ) exit
!       end do
!       if ( j >= i ) then
!         call rfftw_f77_destroy_plan ( antennaPatterns(i)%FFTW_forward_plan )
!         call rfftw_f77_destroy_plan ( antennaPatterns(i)%FFTW_backward_plan )
!       end if
!     end do

    ! Destroy the rest of the database
    do i = 1, size(AntennaPatterns)
      do j = 1, size(AntennaPatterns(i)%signals)
        if (associated(AntennaPatterns(i)%signals(j)%channels)) then
          call deallocate_test ( AntennaPatterns(i)%signals(j)%channels, &
            & "AntennaPatterns(?)%signals(?)%channels", ModuleName )
        end if
      end do
      s = size(AntennaPatterns(i)%signals) * &
        & storage_size(AntennaPatterns(i)%signals) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(AntennaPatterns(i)%signals(1)), addr)
      deallocate ( AntennaPatterns(i)%signals, STAT=status )
      call test_deallocate ( status, moduleName, 'AntennaPatterns(?)%signals', &
        & s, address=addr )
      call deallocate_test ( AntennaPatterns(i)%aaap, &
        & "AntennaPatterns(?)%aaap", moduleName )
      call deallocate_test ( AntennaPatterns(i)%d1aap, &
        & "AntennaPatterns(?)%D1aap", moduleName )
      call deallocate_test ( AntennaPatterns(i)%d2aap, &
        & "AntennaPatterns(?)%D2aap", moduleName )
    end do ! i
    s = size(AntennaPatterns) * storage_size(AntennaPatterns) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(AntennaPatterns(1)), addr)
    deallocate ( AntennaPatterns, stat=status )
    call test_deallocate ( status, moduleName, 'AntennaPatterns', s, address=addr )

  end subroutine Destroy_Ant_Patterns_Database

  ! --------------------------------  Dump_Antenna_Patterns_Database  -----
  subroutine Dump_Antenna_Patterns_Database ( where )
    use Dump_0, only: Dump
    use MoreTree, only: StartErrorMessage
    use Output_m, only: Blanks, Output

    integer, intent(in), optional :: Where   ! Tree node index

    integer :: I, J                ! Subscripts, loop inductors
    if ( associated(antennaPatterns) ) then
      call output ( 'Antenna Patterns: SIZE = ' )
      call output ( size(antennaPatterns), advance='yes' )
      do i = 1, size(AntennaPatterns)
        call output ( i, 4 )
        call output ( ':    Signal = ' )
        do j = 1, size(antennaPatterns(i)%signals)
          if ( j > 1 ) &
          call blanks ( 18 )
          call displaySignalName ( antennaPatterns(i)%signals(j), advance='yes' )
        end do
        call output ( ' Lambda = ' )
        call output ( antennaPatterns(i)%lambda, advance='yes' )
        call dump ( antennaPatterns(i)%aaap, name='Aaap' )
        call dump ( antennaPatterns(i)%d1aap, name='D1aap' )
        call dump ( antennaPatterns(i)%d2aap, name='D2aap' )
      end do ! i
    else
      if ( present(where) ) call startErrorMessage ( where )
      call output ( 'No antenna patterns database to dump.', advance='yes' )
    end if
  end subroutine Dump_Antenna_Patterns_Database

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module AntennaPatterns_m

! $Log$
! Revision 2.18  2015/03/28 01:58:02  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.17  2014/09/05 18:38:01  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.16  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.15  2011/05/09 17:43:04  pwagner
! Converted to using switchDetail
!
! Revision 2.14  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.13  2009/05/13 20:03:01  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.12  2007/10/03 23:58:26  vsnyder
! Add 'where' for tracing
!
! Revision 2.11  2007/06/25 20:33:02  vsnyder
! Add FFTW plans, as comments, in case we ever need them
!
! Revision 2.10  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.9  2004/05/29 02:49:29  vsnyder
! Simplifications from using DisplaySignalName
!
! Revision 2.8  2004/05/26 23:54:14  vsnyder
! Don't dump the database if it's not allocated
!
! Revision 2.7  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.6  2003/05/10 22:20:57  livesey
! Tried to calm down -g1..
!
! Revision 2.5  2003/02/07 01:56:53  vsnyder
! Move a USE down
!
! Revision 2.4  2002/10/08 17:08:01  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.3  2002/09/06 22:31:06  vsnyder
! Cosmetic changes
!
! Revision 2.2  2002/02/22 00:50:08  bill
! forced ntr=4096--wgr
!
! Revision 2.1  2002/02/14 18:39:03  livesey
! Fixed bug with single channel antenna patterns
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.18  2001/05/18 00:00:52  livesey
! Working version.  Still rewinds in the read, but I think
! the file format will have to dictate that for a while.
!
! Revision 1.17  2001/05/17 19:59:38  livesey
! Now pads the arrays to next power of two with zeros after the end
! of the input data.
!
! Revision 1.16  2001/05/16 23:04:20  livesey
! New version, now uses Signal_T, instead of a string
!
! Revision 1.15  2001/05/04 00:48:48  livesey
! Let Destroy database quit if nothing to destroy
!
! Revision 1.14  2001/05/03 22:02:47  vsnyder
! Insert copyright notice and some comments
!
! Revision 1.13  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.12  2001/04/25 23:52:59  livesey
! Added implicit none
!
! Revision 1.11  2001/04/21 01:21:11  vsnyder
! Fix a memory leak
!
! Revision 1.10  2001/04/09 23:45:03  livesey
! Files now two columns rather than 6
!
! Revision 1.9  2001/04/06 00:49:13  vsnyder
! Remove unused variable declarations
!
! Revision 1.8  2001/04/06 00:48:30  vsnyder
! And yet another typo
!
! Revision 1.7  2001/04/06 00:41:50  vsnyder
! Fix another typo
!
! Revision 1.6  2001/04/06 00:21:32  vsnyder
! Fix a typo
!
! Revision 1.5  2001/04/05 22:54:12  vsnyder
! Make array components of AntennaPatterns_T 2-D
!
! Revision 1.4  2001/04/05 00:07:57  vsnyder
! Correct spelling of 'Antenna Pattern'
!
! Revision 1.3  2001/03/30 23:19:06  vsnyder
! Shorten overly-long (standard-violating) subroutine name
!
! Revision 1.2  2001/03/30 20:35:23  zvi
! *** empty log message ***
!
! Revision 1.1  2001/03/30 02:37:44  vsnyder
! Initial Commit
!
