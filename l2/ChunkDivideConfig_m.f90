! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ChunkDivideConfig_m

  ! This module provides ChunkDivideConfig_t, ChunkDivideConfig, and
  ! its dumpers.
  ! The main purpose is to avoid a circular dependence between
  ! ChunkDivide_m and Dump_Command.

  use INTRINSIC, only: L_NONE, PHYQ_INVALID
  use MLSKinds, only: RP
  use MLSSTRINGLISTS, only: SWITCHDETAIL
  use TOGGLES, only: SWITCHES

  implicit NONE
  private

  public :: ChunkDivideConfig_t, Dump, Dump_Config, Dump_criticalSignals

  ! This type is filled by the l2cf and describes the configuration of the
  ! chunk division process.
  type ChunkDivideConfig_T
    integer :: method = l_none          ! See below.
    real(rp) :: maxLength = 0           ! Maximum length of chunks
    integer :: maxLengthFamily = PHYQ_Invalid ! PHYQ_Angle etc.
    integer :: noChunks = 0             ! Number of chunks for [fixed]
    real(rp) :: overlap = 0.0           ! Desired length of overlaps
    real(rp) :: lowerOverlap = 0.0      ! Desired length of lower overlaps
    real(rp) :: upperOverlap = 0.0      ! Desired length of lower overlaps
    integer :: overlapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
    integer :: lowerOverlapFamily = PHYQ_Invalid
    integer :: upperOverlapFamily = PHYQ_Invalid
    integer :: noSlaves = 0             ! Number of slave nodes [even]
    integer :: homeModule = l_none      ! Which module to consider [orbital]
    real(rp) :: homeGeodAngle = 0.0     ! Aim for one chunk to start here [orbital]
    logical   :: scanLLSet = .false.    ! True if scan lower limit should be used
    logical   :: scanULSet = .false.    ! True if scan upper limit should be used
    real(rp), dimension(2) :: scanLowerLimit ! Range for bottom of scan
    real(rp), dimension(2) :: scanUpperLimit ! Range for top of scan
    real(rp) :: maxOrbY = -1.0                 ! Maximum out of plane distance allowed <=0.0 default
    character(len=128) :: criticalBands = ' ' ! Which bands must be scanning
    integer   :: criticalModules = l_none ! Which modules must be scanning
    logical   :: chooseCriticalSignals = .true. ! Use criticalModules, Bands?
    character(len=160), dimension(:), pointer &
      & :: criticalSignals => null()    ! Which signals must be on
    real(rp)  :: maxGap = 0.0           ! Length of time/MAFs/orbits allowed for gap
    integer   :: maxGapFamily = PHYQ_Invalid ! PHYQ_MAF, PHYQ_Time etc.
    logical   :: skipL1BCheck = .false. ! Don't check for l1b data probs
    logical   :: crashIfPhiNotMono = .false. ! If l1b contains non-monotonic phi
    logical   :: allowPriorOverlaps = .true. ! Use MAFs before start time
    logical   :: allowPostOverlaps = .true. ! Use MAFs after end time
    logical   :: saveObstructions = .true. ! Save obstructions for Output_Close
    logical   :: DACSDeconvolved = .true. ! Don't need to do this in level 2
    integer   :: numPriorOverlaps = 0   ! How many profiles before processingRanges
    integer   :: numPostOverlaps = 0    ! How many profiles after processingRanges
    integer   :: Where = 0              ! in the l2cf tree it was defined
  end type ChunkDivideConfig_T

  type(ChunkDivideConfig_T), public, save :: CHUNKDIVIDECONFIG

  interface dump
    module procedure DUMP_CONFIG
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

contains ! ===================================  Public Procedures  =====

  ! -------------------------------------------- Dump_Config -----
  subroutine Dump_Config ( config )

    use Intrinsic, only: Lit_Indices, PHYQ_Indices
    use LEXER_CORE, only: PRINT_SOURCE
    use Output_m, only: Output
    use String_Table, only: Display_String
    use TREE, only: WHERE

    ! Args
    type(ChunkDivideConfig_T), intent(in) :: Config

    ! Executable code
    call print_source ( where(config%where), &
      & before='ChunkDivide configuration at ', advance='yes' )
    call display_string ( lit_indices(Config%method), &
      &             strip=.true., before='  method ', advance='yes' )
    call output ( config%maxLength, before='  max Length ', advance='yes' )
    call display_string ( phyq_indices(Config%maxLengthFamily), &
      &             strip=.true., before='  max Length Family ', advance='yes' )
    call output ( config%noChunks, before='  num chunks ', advance='yes' )
    call output ( config%overlap, before='  overlap ', advance='yes' )
    call display_string ( phyq_indices(Config%overlapFamily), &
      &             strip=.true., before='  overlap Family ', advance='yes' )
    call output ( config%loweroverlap, before='  lower overlap ', advance='yes' )
    call display_string ( phyq_indices(Config%loweroverlapFamily), &
      &             strip=.true., before='  lower overlap Family ', advance='yes' )
    call output ( config%upperoverlap, before='  upper overlap ', advance='yes' )
    call display_string ( phyq_indices(Config%upperoverlapFamily), &
      &             strip=.true., before='  upper overlap Family ', advance='yes' )
    call output ( config%noSlaves, before='  num slaves ', advance='yes' )
    call display_string ( lit_indices(Config%homeModule), &
      &             strip=.true., before='  home module ', advance='yes' )
    call output ( config%homeGeodAngle, before='  home Geod Angle ', advance='yes' )
    call output ( config%scanLLSet, before='  set scan lower limit? ', advance='yes' )
    if ( config%scanLLSet ) then
      call output ( '  Bottom scan range ' )
      call output ( config%scanLowerLimit, advance='yes' )
    end if
    call output ( config%scanULSet, before='  set scan upper limit? ', advance='yes' )
    if ( config%scanULSet ) then
      call output ( '  Top scan range ' )
      call output ( config%scanUpperLimit, advance='yes' )
    end if
    call output ( config%maxOrbY, before='  max out-of-plane distance ', &
      & advance='yes' )
    call display_string ( lit_indices(Config%criticalModules), &
      &             strip=.true., before='  critical modules ', advance='yes' )
    call output ( '  critical bands ' )
    call output ( trim(config%criticalBands), advance='yes' )
    call output ( config%chooseCriticalSignals, &
      & before='  use critical modules to choose critical signals? ', advance='yes' )
    call output ( config%maxGap, before='  max gap ', advance='yes' )
    call display_string ( phyq_indices(Config%maxGapFamily), &
      &             strip=.true., before='  max Gap Family ', advance='yes' )
    call output ( config%skipL1BCheck, before='  skip check of l1b files ', &
      & advance='yes' )
    call output ( config%allowPriorOverlaps, &
      & before='  allow overlaps to prior day? ',advance='yes' )
    call output ( config%allowPostOverlaps, &
      & before='  allow overlaps to next day? ', advance='yes' )
    call output ( config%saveObstructions, before='  save obstructions? ', &
      & advance='yes' )
    call output ( config%DACSDeconvolved, before='  DACS already deconvolved? ', &
      & advance='yes' )
    call Dump_criticalSignals(config%criticalSignals)
  end subroutine Dump_Config

  ! -------------------------------------------- Dump_CriticalSignals -----
  subroutine Dump_CriticalSignals(criticalSignals)

    use Output_M, only: OUTPUT

    character(len=160), dimension(:), pointer &
      & :: criticalSignals       ! Which signals must be on

    ! Local variables
    integer :: i                        ! Loop counter

    ! Executable code
    if ( associated ( criticalSignals ) ) then
      if ( size(criticalSignals) == 0 ) then
        call output ( 'criticalSignals is a zero size array.', advance='yes' )
      else
        call output ( 'Dumping ' )
        call output ( size(criticalSignals) )
        call output ( ' criticalSignals:', advance='yes' )
        do i = 1, size(criticalSignals)
          call output ( i )
          call output ( ': ' )
          call output ( trim(criticalSignals(i)), advance='yes' )
        end do
      end if
    else
      call output ( 'critical Signals is not associated.', advance='yes')
    end if
  end subroutine Dump_CriticalSignals

! ===========================================  Private Procedures  =====

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ChunkDivideConfig_m

! $Log$
! Revision 2.3  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.2  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.1  2013/09/21 00:26:14  vsnyder
! Initial commit; avoid circular dependence between ChunkDivide and DumpCommand
!
