! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Timings              !  Timings for the MLSL2 program sections
!=============================================================================

  use INTRINSIC, only: L_HOURS, L_MINUTES, L_SECONDS
  use L2PARINFO, only: PARALLEL
  USE MLSL2Options, only: SECTIONTIMINGUNITS
  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
  USE MLSStrings, only: GetStringElement, LowerCase, StringElementNum
  use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
  use Time_M, only: Time_Now, Use_Wall_Clock

  IMPLICIT NONE

  public :: SECTION_TIMES, TOTAL_TIMES, &
    & add_to_directwrite_timing, add_to_phase_timing, &
    & add_to_retrieval_timing, add_to_section_timing, &
    & dump_section_timings, run_start_time
  private

  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! This module simply contains initial settings and accumulated values.

  ! The following public settings are stored here; they may be set by MLSL2
  logical          :: SECTION_TIMES = .false.  ! Show times in each section
  logical          :: TOTAL_TIMES = .false.    ! Show total times from start

  logical, private   :: COUNTEMPTY = .false.     ! Any sections named ' '?
  logical, private   :: ALLOWUNKNOWNNAMES = .false.  ! Any unknown section names
  integer, parameter :: PHASENAMESLEN = 400
  integer, parameter :: MAXNUMPHASES = 40

  character*(*), parameter           :: section_names = &
    & 'main,open_init,global_settings,signals,spectroscopy,' // &
    & 'read_apriori,chunk_divide,construct,fill,retrieve,join,' // &
    & 'directwrite,algebra,output'
  ! This should be the number of elements in the above ---------------
  integer, parameter                 :: num_section_times = 14    ! <--|

  character*(*), parameter           :: retrieval_names = &
    & 'newton_solver,cholesky_factor,cholesky_solver,cholesky_invert,' // &
    & 'baseline,hybrid,polar_linear,switching_mirror,' // &
    & 'full_fwm,fullcloud_fwm,scan_fwm,twod_scan_fwm,linear_fwm,' // &
    & 'low_cloud,high_cloud,sids,form_normeq,tikh_reg'
  ! This should be the number of elements in the above -----------------
  integer, parameter                 :: num_retrieval_times = 18  ! <--|

  character*(*), parameter           :: directwrite_names = &
    & 'writing,waiting'
  ! This should be the number of elements in the above ----------------
  integer, parameter                 :: num_directwrite_times = 2  ! <--|
  ! dimension of the following is +2 to allow possible unknown
  ! section names and unknown retrieval names
  real, dimension(num_section_times+num_retrieval_times+num_directwrite_times+2), &
    & save                           :: section_timings = 0.
  ! integer, parameter                 :: unknown_section = &
  !   &                               num_section_times+num_retrieval_times + 1
  ! integer, parameter                 :: unknown_retrieval = unknown_section + 1
  character(len=PHASENAMESLEN), save :: phaseNames = ' '
  integer, save :: num_phases = 0
  real, dimension(MAXNUMPHASES), save :: phase_timings = 0.
  real, save :: run_start_time = 0.

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------  add_to_directwrite_timing  -----
  subroutine add_to_directwrite_timing( section_name, t1 )
  ! Add current elapsed directwrite section time to total so far for section_name

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the dw. sect_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    integer                     :: elem_offset
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
      if ( present(t1) ) myLastTime = t1
      elem_offset = num_section_times + num_retrieval_times

      elem = StringElementNum(directwrite_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_directwrite_times ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find directwrite section name ' // section_name // &
        & ' among list ' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(elem_offset+elem) = &
          & section_timings(elem_offset+elem) + t2 - myLastTime
      endif
      myLastTime = t2
      if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_directwrite_timing

  ! -----------------------------------------  add_to_phase_timing  -----
  subroutine add_to_phase_timing( phase_name, t1 )
  ! Add current elapsed phase time to total so far for last phase

  ! Formal arguments
    character(LEN=*), intent(in):: phase_name   ! One of the phases, e.g. 'core'
    real, optional, intent(inout)  :: t1        ! Prior time_now, then current

  ! Private
    integer                     :: elem
    integer, save               :: myLastElem = 0
    real                        :: t2
    real, save                  :: myLastTime = 0.

    ! Executable
      if ( present(t1) ) myLastTime = t1
      elem = StringElementNum(LowerCase(trim(phaseNames)), &
        & LowerCase(phase_name), countEmpty)
      if ( elem < 1 .and. phase_name /= ' ' ) then
        num_phases = num_phases + 1
        if ( num_phases > MAXNUMPHASES ) &
          & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Too many phases--unable add to phase name ' // phase_name // &
          & ' to list ' // trim(phaseNames) )
        phaseNames = trim(phaseNames) // ',' // trim(phase_name)
        elem = num_phases
      endif
      call time_now ( t2 )
      if ( myLastElem > 0 ) then
        phase_timings(myLastElem) = &
          & phase_timings(myLastElem) + t2 - myLastTime
      endif
      myLastTime = t2
      myLastElem = elem
      if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_phase_timing

  ! -----------------------------------------  add_to_retrieval_timing  -----
  subroutine add_to_retrieval_timing( section_name, t1 )
  ! Add current elapsed retrieval section time to total so far for section_name

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the retr. sect_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
      if ( present(t1) ) myLastTime = t1

      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find retrieval section name ' // section_name // &
        & ' among list ' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(num_section_times+elem) = &
          & section_timings(num_section_times+elem) + t2 - myLastTime
      endif
      myLastTime = t2
      if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_retrieval_timing

  ! -----------------------------------------------  add_to_section_timing  -----
  subroutine add_to_section_timing( section_name, t1 )
  ! Add current elapsed section time to total so far for section_name
  ! (or possibly, one of the retrieval sections)

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the section_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
    if ( present(t1) ) myLastTime = t1
    elem = StringElementNum(section_names, LowerCase(section_name), countEmpty)
    if ( elem < 1 .or. elem > num_section_times ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Unable to find section name ' // section_name // &
      & ' among list ' // section_names // ',' // retrieval_names )
    else
      call time_now ( t2 )
      section_timings(elem) = section_timings(elem) + t2 - myLastTime
    endif
    myLastTime = t2
    if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_section_timing

  ! -----------------------------------------------  dump_section_timings  -----
  subroutine dump_section_timings
  ! (1)dump accumulated elapsed timings for section_names
  ! (2) breakdown timings by subsection in retrieve, directwrite
  ! (3) breakdown timings by phase
  ! Private
    character(LEN=*), parameter     :: TIMEFORMSMALL = '(F10.2)'
    character(LEN=*), parameter     :: TIMEFORMBIG = '(1PE10.2)'
    character(LEN=*), parameter     :: PCTFORM='(F10.0)'
    logical, parameter              :: PRINTCHUNKNUMWITHPHASES=.true.
    integer                         :: elem
    integer                         :: joinElem
    integer                         :: dwElem
    integer                         :: retrElem
    character(LEN=16)               :: section_name   ! One of the section_names
    real                            :: final
    real                            :: total
    real                            :: retrFinal
    real                            :: retrTotal
    real                            :: dwFinal
    real                            :: dwTotal
    real                            :: phaseTotal
    real                            :: percent
    real                            :: elem_time
    character(LEN=LEN(TIMEFORMBIG)) :: TIMEFORM
    logical                         :: Unknown_nonzero
    integer                         :: num_elems
    integer                         :: timeDivisor

  ! Executable
    select case ( sectionTimingUnits )
    case ( l_seconds )
      timeDivisor = 1
    case ( l_minutes )
      timeDivisor = 60
    case ( l_hours )
      timeDivisor = 3600
    case default
      call output ( 'Unrecognized section timing unit : ', advance='no' )
      call output ( sectionTimingUnits, advance='yes' )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Should be among list l_(seconds, minutes, hours)' )
    end select
    Unknown_nonzero = .false.  ! (section_timings(unknown_section) > 0.)
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Level 2 section timings : ', advance='yes' )
    call blanks ( 8, advance='no' )
    call printTaskType
    call output ( '==========================================', advance='yes' )
    call output ( 'section name ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'percent ', advance='yes' )

    ! A trick:
    ! Adjust Join section timing to exclude directwrite timings
    ! (otherwise they would be counted twice)
    joinElem = StringElementNum(section_names, 'join', countEmpty)
    dwElem = StringElementNum(section_names, 'directwrite', countEmpty)
    if ( joinElem > 0 .and. dwElem > 0 ) section_timings(joinElem) = &
      &         section_timings(joinElem) - section_timings(dwElem)

    total = sum(section_timings(1:num_section_times)) ! + &
     ! & section_timings(unknown_section)

    ! Another trick:
    ! The DirectWrite section doesn't automatically include the waiting time
    ! (due to lazy coding in Join) so add it in now
    elem = StringElementNum(directwrite_names, 'waiting', countEmpty)
    section_timings(dwElem) = section_timings(dwElem) + &
      & section_timings(num_section_times+num_retrieval_times+elem)
    call time_now ( final )
    final = max(final-run_start_time, total)
    if ( final/timeDivisor <= 0.0 ) final = 1.0       ! Just so we don't divide by 0
    if ( final/timeDivisor < 0.5 .or. final/timeDivisor > 99999.99 ) then
      TIMEFORM = TIMEFORMBIG
    else
      TIMEFORM = TIMEFORMSMALL
    endif
    ! if ( Unknown_nonzero ) then
    !  num_elems = num_section_times + 1
    ! else
    !  num_elems = num_section_times
    ! endif
    do elem = 1, num_section_times
      ! if ( elem > num_section_times ) then
      !  elem_time = section_timings(unknown_section)
      !  section_name = 'unknown name'
      ! else
        elem_time = section_timings(elem)
        call GetStringElement(section_names, section_name, elem, countEmpty)
      ! endif
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    enddo
    call output ( '==========================================', advance='yes' )
    percent = 100 * total / final
    call blanks ( 2, advance='no' )
    call output ( '(subtotal)', advance='no' )
    call blanks ( 5, advance='no' )
    call output ( total/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    percent = 100 * (final-total) / final
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( (final-total)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    call output ( '==========================================', advance='yes' )
    percent = 100
    call blanks ( 4, advance='no' )
    call output ( '(final)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( final/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='no' )
    call blanks ( 4, advance='no' )
    call printTaskType
    ! Subdivision of Retrieval section
    retrElem = StringElementNum(section_names, 'retrieve', countEmpty)
    if ( retrElem == 0 ) then
      call output ( '(Illegal section name--spelling?) ', advance='yes' )
      return
    endif
    retrFinal = section_timings(retrElem) 
    if ( retrFinal == 0.0 ) then
      call output ( '(Retrieval section number ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( retrFinal, advance='yes' )
      call output ( '(No Retrieval section timings breakdown) ', advance='yes' )
    else
      call output ( '==========================================', advance='yes' )
      call blanks ( 8, advance='no' )
      call output ( 'Retrieval section timings : ', advance='yes' )
      call output ( '==========================================', advance='yes' )
      call output ( 'subsection name ', advance='no' )
      call blanks ( 8, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'percent of total retrieval time', advance='yes' )
      retrTotal = sum(section_timings(1+num_section_times:num_section_times+num_retrieval_times)) ! - &
        ! & section_timings(unknown_section)
      Unknown_nonzero = .false.  ! (section_timings(unknown_retrieval) > 0.)
      ! if ( Unknown_nonzero ) then
      !  num_elems = num_retrieval_times + 1
      !  else
      ! num_elems = num_retrieval_times
      ! endif
      do elem = 1, num_retrieval_times  ! num_elems
        ! if ( elem > num_retrieval_times ) then
        !  elem_time = section_timings(unknown_retrieval)
        !  section_name = 'unknown name'
        !else
          elem_time = section_timings(num_section_times+elem)
          call GetStringElement(retrieval_names, section_name, elem, countEmpty)
        !endif
        percent = 100 * elem_time / retrFinal
        call output ( section_name, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
      enddo
      call blanks ( 3, advance='no' )
      call output ( '(others)', advance='no' )
      call blanks ( 6, advance='no' )
      call output ( (retrFinal-retrTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 100*(retrFinal-retrTotal)/retrFinal, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    endif

    ! Subdivision of DirectWrite section
    elem = StringElementNum(section_names, 'directwrite', countEmpty)
    if ( elem == 0 ) then
      call output ( '(Illegal section name--spelling?) ', advance='yes' )
      return
    endif
    dwFinal = section_timings(elem) 
    if ( dwFinal == 0.0 ) then
      call output ( '(DirectWrite section number ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( dwFinal, advance='yes' )
      call output ( '(No DirectWrite section timings breakdown) ', &
        & advance='yes' )
    else
      call output ( '==========================================', advance='yes' )
      call blanks ( 8, advance='no' )
      call output ( 'DirectWrite section timings : ', advance='yes' )
      call output ( '==========================================', advance='yes' )
      call output ( 'subsection name ', advance='no' )
      call blanks ( 8, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'percent of total DirectWrite time', advance='yes' )
      dwTotal = sum(section_timings(1+num_section_times+num_retrieval_times:)) 
      do elem = 1, num_directwrite_times
        elem_time = section_timings(num_section_times+num_retrieval_times+elem)
        call GetStringElement(directwrite_names, section_name, elem, countEmpty)
        percent = 100 * elem_time / dwFinal
        call output ( section_name, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, &
          & advance='no' )
        call blanks ( 2, advance='no' )
        call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
      enddo
      call blanks ( 3, advance='no' )
      call output ( '(others)', advance='no' )
      call blanks ( 6, advance='no' )
      call output ( (dwFinal-dwTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, &
        & advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 100*(dwFinal-dwTotal)/dwFinal, FORMAT=PCTFORM, &
        & LOGFORMAT=PCTFORM, advance='yes' )
    endif

    ! Subdivision of Phases
    if ( num_phases == 0 ) then
      call output ( '(No phase timings breakdown) ', advance='yes' )
      return
    endif
    ! A trick! Add final elapsed time to last phase
    call add_to_phase_timing(' ')
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Phase timings : ', advance='yes' )
    call output ( '==========================================', advance='yes' )
    call output ( 'phase name ', advance='no' )
    call blanks ( 11, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 5, advance='no' )
    call output ( 'percent of total time', advance='yes' )
    phaseTotal = 0.
    do elem = 1, num_phases
      elem_time = phase_timings(elem)
      call GetStringElement(trim(phaseNames), section_name, elem, countEmpty)
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='no' )
      if ( PRINTCHUNKNUMWITHPHASES ) then
        call printTaskType
      else
        call output ( ' ', advance='yes' )
      endif
      phaseTotal = phaseTotal + elem_time
    enddo
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( (final-phaseTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( 100*(final-phaseTotal)/final, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )

    contains
    ! Internal subprograms
    subroutine printTaskType
      if ( parallel%master ) then
        call output ( '(Master Task) ', advance='yes' )
      elseif ( parallel%slave .and. .not. parallel%fwmParallel ) then
        call output ( '(Slave: chunk ', advance='no' )
        call output ( parallel%ChunkNo, advance='no' )
        call output ( ' ) ', advance='yes' )
      else
        call output ( '(Serial Task) ', advance='yes' )
      endif
    end subroutine printTaskType
  end subroutine dump_section_timings

!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END MODULE MLSL2Timings
!=============================================================================

!
! $Log$
! Revision 2.21  2004/02/10 19:26:23  pwagner
! Fixed bug in phase names
!
! Revision 2.20  2004/01/14 18:49:58  vsnyder
! Stuff to support the Algebra section
!
! Revision 2.19  2003/12/05 00:50:44  pwagner
! Numerous small bugixes related to wall clock time use
!
! Revision 2.18  2003/10/23 22:20:16  pwagner
! A few bugfixes in dump_section_timings
!
! Revision 2.17  2003/10/22 21:19:14  pwagner
! Added timings breakdown by phases
!
! Revision 2.16  2003/10/20 18:21:45  pwagner
! Timings breakdown added for directWrite
!
! Revision 2.15  2003/08/11 23:24:48  pwagner
! Chunk no. printed if slave task
!
! Revision 2.14  2003/06/09 22:51:36  pwagner
! Renamed scan_divide to chunk_divide in timings table
!
! Revision 2.13  2003/02/27 21:56:07  pwagner
! Passes LOGFORMAT along with FORMAT
!
! Revision 2.12  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2002/09/24 18:15:12  pwagner
! Prepared to allow unknown names; add_to_section_timing now like retrieval
!
! Revision 2.10  2002/09/19 19:07:05  vsnyder
! Only update t1 if it's present!
!
! Revision 2.9  2002/09/18 23:56:01  vsnyder
! Call time_now at end of add_to_retrieval_timing
!
! Revision 2.8  2002/07/23 00:06:05  pwagner
! No upper-case allowed in section names
!
! Revision 2.7  2002/07/22 22:53:10  pwagner
! Added names of 2d scan model, form norm eq, and tikh reg to retrieval
!
! Revision 2.6  2001/11/27 23:34:49  pwagner
! Split forward model timings into four types
!
! Revision 2.5  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.4  2001/10/01 23:30:50  pwagner
! Fixed bug in spelling cholesky_solver
!
! Revision 2.3  2001/10/01 22:54:22  pwagner
! Added subsection timings for Retrieval section
!
! Revision 2.2  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.1  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
