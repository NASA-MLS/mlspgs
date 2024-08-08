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

  use Intrinsic, only: L_None, Phyq_Invalid
  use MLSKinds, only: Rp

  implicit none
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
    integer :: module = l_none  
    real(rp) :: homeGeodAngle = 0.0     ! Aim for one chunk to start here [orbital]
    logical   :: scanLLSet = .false.    ! True if scan lower limit should be used
    logical   :: scanULSet = .false.    ! True if scan upper limit should be used
    real(rp), dimension(2) :: scanLowerLimit ! Range for bottom of scan
    real(rp), dimension(2) :: scanUpperLimit ! Range for top of scan
    real(rp) :: maxOrbY = -1.0                 ! Maximum out of plane distance allowed <=0.0 default
    real(rp) :: maxScVel = 1000000000.         ! Maximum s/c velocity (abs val)
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

  type(ChunkDivideConfig_T), public, save :: ChunkDivideConfig

  interface Dump
    module procedure Dump_config
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

contains ! ===================================  Public Procedures  =====

  ! -------------------------------------------- Dump_Config -----
  subroutine Dump_Config ( config )

    use HighOutput, only: AddRow, AddRow_Divider, AddRow_Header, &
      & OutputTable, StartTable
    use Intrinsic, only: Lit_Indices, PHYQ_Indices
    use Lexer_core, only: Print_Source
    use MLSStrings, only: LowerCase
    use String_Table, only: Get_String
    use Tree, only: Where

    ! Args
    type(ChunkDivideConfig_T), intent(in) :: Config

    ! Executable code
    call print_source ( where(config%where), &
      & before='ChunkDivide configuration at ', advance='yes' )
    call startTable
    call addRow_header ( 'ChunkDivide configuration', 'c' )
    call addRow_divider ( '-' )
    call addRow ( 'Method                    ', trim(string_function(Config%method, 'lit')  ) )
    call addRow ( 'Max length                ', config%maxLength )
    call addRow ( 'Max length   family       ', trim(string_function(Config%maxLengthFamily, 'phy')  ) )
    
    call addRow ( 'Num chunks                ', config%noChunks )
    call addRow ( 'Overlap                   ', config%overlap )
    call addRow ( 'Overlap  family           ', trim(string_function(Config%overlapFamily, 'phy')  ) )
    call addRow ( 'Lower overlap             ', config%loweroverlap )
    call addRow ( 'Its  family               ', trim(string_function(Config%loweroverlapFamily, 'phy')  ) )
    call addRow ( 'Upper overlap             ', config%upperoverlap )
    call addRow ( 'Its  family               ', trim(string_function(Config%upperoverlapFamily, 'phy')  ) )
    call addRow ( 'Num slaves                ', config%noSlaves )
    call addRow ( 'Home Module               ', trim(string_function(Config%homeModule, 'lit')  ) )
    call addRow ( 'Home Geod Ang             ', config%homeGeodAngle )
    call addRow ( 'Set Scan Lower Limit?     ', config%ScanLLSet )
    if ( config%scanLLSet ) &
  & call addRow ( 'Bottom Scan Range         ', config%scanlowerLimit )

    call addRow ( 'Set Scan Upper Limit?     ', config%ScanULSet )
    if ( config%scanULSet ) &
  & call addRow ( 'Top Scan Range            ', config%scanUpperLimit )
    call addRow ( 'Max Out-of-plane Dist     ', config%MaxOrbY )
    call addRow ( 'Max s/c vel               ', config%MaxScVel )
    call addRow ( 'Critical Modules          ', trim(string_function(Config%criticalModules, 'lit')  ) )
    call addRow ( 'Critical Bands            ', trim(Config%criticalBands) )
    call addRow ( 'Use Crit. Modules?        ', config%chooseCriticalSignals )
    call addRow ( 'Max Gap                   ', config%MaxGap )
    call addRow ( 'Max Gap  family           ', trim(string_function(Config%maxGapFamily, 'phy')  ) )
    call addRow ( 'Skip L1B Check?           ', config%skipL1BCheck )
    call addRow ( 'Allow Prior Overlaps?     ', config%allowPriorOverlaps )
    call addRow ( 'Allow Next Day Overlaps?  ', config%allowPostOverlaps )
    call addRow ( 'Save Obstructions?        ', config%saveObstructions )
    call addRow ( 'DACS Already Deconvolved? ', config%DACSDeconvolved )
    call outputTable ( sep='|', border='-' )
    ! Critical signals?
    call Dump_criticalSignals(config%criticalSignals)
  contains
    function string_function ( arg, typ ) result ( the_string )
      ! Returns the value get_string computes
      ! Args
      integer, intent(in)           :: arg
      character(len=*), intent(in)  :: typ
      character(len=128)            :: the_string
      ! Executable
      select case (lowercase(typ(1:3)))
      case ( 'lit' )
        call Get_string ( lit_indices(arg), the_string, strip=.true. )
      case ( 'phy' )
        call Get_string ( phyq_indices(arg), the_string, strip=.true. )
      case default
        call Get_string ( arg, the_string, strip=.true. )
      end select
    end function string_function
  end subroutine Dump_Config

  ! -------------------------------------------- Dump_CriticalSignals -----
  subroutine Dump_CriticalSignals(criticalSignals)
    ! Some day we'll use the Table apis from highOutput (see Dump_Config)

    use Output_M, only: Output

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
! Revision 2.7  2024/08/08 20:43:39  pwagner
! s/c velocity above maxScVel means an invalid MAF now
!
! Revision 2.6  2019/07/09 20:46:35  pwagner
! Use Table ccells to Dump_Config
!
! Revision 2.5  2017/09/14 23:17:26  pwagner
! Added module component to ChunkDivideConfig_T
!
! Revision 2.4  2014/08/06 23:26:45  vsnyder
! Remove USE for Switches and SwitchDetail, which are not referenced
!
! Revision 2.3  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.2  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.1  2013/09/21 00:26:14  vsnyder
! Initial commit; avoid circular dependence between ChunkDivide and DumpCommand
!
