! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Timings              !  Timings for the MLSL2 program sections
!=============================================================================

  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
  USE MLSStrings, only: GetStringElement, LowerCase, StringElementNum
  use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module simply contains initial settings and accumulated values.

  logical            :: SECTION_TIMES=.false.  ! Show times in each section
  logical            :: TOTAL_TIMES=.false.    ! Show total times from start
  logical, private   :: COUNTEMPTY=.false.     ! Any sections named ' '?

  character*(*), parameter           :: section_names = &
    & 'open_init, global_settings, signals, spectroscopy, read_apriori, ' // &
    & 'scan_divide, construct, fill, retrieve, join, output'
  integer, parameter                 :: num_section_times = 11
  real, dimension(num_section_times) :: section_timings = 0.

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  add_to_section_timing  -----
  subroutine add_to_section_timing( section_name, t1 )
  ! Add current elapsed section time to total so far for section_name

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the section_names
    real, intent(in)            :: t1             ! Prior cpu_time 

  ! Private
    integer                     :: elem
    real                        :: t2

  ! Executable
    elem = StringElementNum(section_names, LowerCase(section_name), countEmpty)
    if ( elem < 1 .or. elem > num_section_times ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find section name ' // section_name // &
        & 'among list ' // section_names)
    else
      call cpu_time ( t2 )
      section_timings(elem) = section_timings(elem) + t2 - t1
    endif
  end subroutine add_to_section_timing

  ! -----------------------------------------------  dump_section_timings  -----
  subroutine dump_section_timings
  ! dump accumulated elapsed timings forsection_names

  ! Private
    integer                     :: elem
    character(LEN=16)           :: section_name   ! One of the section_names
    real                        :: final
    real                        :: total
    real                        :: percent

  ! Executable
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Level 2 section timings : ', advance='yes' )
    call output ( '==========================================', advance='yes' )
    call output ( 'section name ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'percent ', advance='yes' )
    total = sum(section_timings)
    call cpu_time ( final )
    final = max(final, total)
    if ( final <= 0.0 ) final = 1.0       ! Just so we don't divide by 0
    do elem = 1, num_section_times
      call GetStringElement(section_names, section_name, elem, countEmpty)
      percent = 100 * section_timings(elem) / final
      call output ( section_name, advance='no' )
      call blanks ( 8, advance='no' )
      call output ( section_timings(elem), advance='no' )
      call blanks ( 8, advance='no' )
      call output ( percent, advance='yes' )
    enddo
    percent = 100 * (final-total) / final
    call output ( '(others)', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( final-total, advance='no' )
    call blanks ( 8, advance='no' )
    call output ( percent, advance='yes' )
  end subroutine dump_section_timings

!=============================================================================
END MODULE MLSL2Timings
!=============================================================================

!
! $Log$
! Revision 2.1  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
