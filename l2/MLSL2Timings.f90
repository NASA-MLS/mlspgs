! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Timings              !  Timings for the MLSL2 program sections
!=============================================================================

  USE MLSMessageModule, only: MLSMessage, MLSMSG_Error
  USE MLSStrings, only: GetStringElement, LowerCase, StringElementNum
  use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
  use Time_M, only: Time_Now

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module simply contains initial settings and accumulated values.

  logical            :: SECTION_TIMES = .false.  ! Show times in each section
  logical            :: TOTAL_TIMES = .false.    ! Show total times from start
  logical, private   :: COUNTEMPTY = .false.     ! Any sections named ' '?

  character*(*), parameter           :: section_names = &
    & 'main,open_init,global_settings,signals,spectroscopy,' // &
    & 'read_apriori,scan_divide,construct,fill,retrieve,join,output'
  integer, parameter                 :: num_section_times = 12
  character*(*), parameter           :: retrieval_names = &
    & 'newton_solver,cholesky_factor,cholesky_solver,cholesky_invert,' // &
    & 'forward_model,low_cloud,sids'
  integer, parameter                 :: num_retrieval_times = 7
  real, dimension(num_section_times+num_retrieval_times), &
    & save                           :: section_timings = 0.

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  add_to_retrieval_timing  -----
  subroutine add_to_retrieval_timing( section_name, t1 )
  ! Add current elapsed retrieval section time to total so far for section_name

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the retrieval section_names
    real, intent(in)            :: t1             ! Prior time_now 

  ! Private
    integer                     :: elem
    real                        :: t2

  ! Executable
      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find retrieval section name ' // section_name // &
        & ' among list ' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(num_section_times+elem) = &
          & section_timings(num_section_times+elem) + t2 - t1
      endif
  end subroutine add_to_retrieval_timing

  ! -----------------------------------------------  add_to_section_timing  -----
  subroutine add_to_section_timing( section_name, t1 )
  ! Add current elapsed section time to total so far for section_name
  ! (or possibly, one of the retrieval sections)

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the section_names
    real, intent(in)            :: t1             ! Prior time_now 

  ! Private
    integer                     :: elem
    real                        :: t2

  ! Executable
    elem = StringElementNum(section_names, LowerCase(section_name), countEmpty)
    if ( elem < 1 .or. elem > num_section_times ) then
      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
        call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find section name ' // section_name // &
        & ' among list ' // section_names // ',' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(num_section_times+elem) = &
          & section_timings(num_section_times+elem) + t2 - t1
      endif
    else
      call time_now ( t2 )
      section_timings(elem) = section_timings(elem) + t2 - t1
    endif
  end subroutine add_to_section_timing

  ! -----------------------------------------------  dump_section_timings  -----
  subroutine dump_section_timings
  ! dump accumulated elapsed timings forsection_names

  ! Private
    character(LEN=*), parameter     :: TIMEFORMSMALL = '(F10.2)'
    character(LEN=*), parameter     :: TIMEFORMBIG = '(1PE10.2)'
    character(LEN=*), parameter     :: PCTFORM='(F10.0)'
    integer                         :: elem
    character(LEN=16)               :: section_name   ! One of the section_names
    real                            :: final
    real                            :: total
    real                            :: percent
    character(LEN=LEN(TIMEFORMBIG)) :: TIMEFORM

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
    total = sum(section_timings(1:num_section_times))
    call time_now ( final )
    final = max(final, total)
    if ( final <= 0.0 ) final = 1.0       ! Just so we don't divide by 0
    if ( final < 0.5 .or. final > 99999.99 ) then
      TIMEFORM = TIMEFORMBIG
    else
      TIMEFORM = TIMEFORMSMALL
    endif
    do elem = 1, num_section_times
      call GetStringElement(section_names, section_name, elem, countEmpty)
      percent = 100 * section_timings(elem) / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( section_timings(elem), FORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( percent, FORMAT=PCTFORM, advance='yes' )
    enddo
    call output ( '==========================================', advance='yes' )
    percent = 100 * total / final
    call blanks ( 4, advance='no' )
    call output ( '(total)', advance='no' )
    call blanks ( 7, advance='no' )
    call output ( total, FORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, advance='yes' )
    percent = 100 * (final-total) / final
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 7, advance='no' )
    call output ( final-total, FORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, advance='yes' )
    call output ( '==========================================', advance='yes' )
    percent = 100
    call blanks ( 4, advance='no' )
    call output ( '(final)', advance='no' )
    call blanks ( 7, advance='no' )
    call output ( final, FORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, advance='yes' )

    ! Subdivision of Retrieval section
    elem = StringElementNum(section_names, 'retrieve', countEmpty)
    if ( final == 0.0 ) then
      call output ( '(Illegal section name--spelling?) ', advance='yes' )
      return
    endif
    final = section_timings(elem) 
    if ( final == 0.0 ) then
      call output ( '(Retrieval section number ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( final, advance='yes' )
      call output ( '(No Retrieval section timings breakdown) ', advance='yes' )
      return
    endif
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Retrieval section timings : ', advance='yes' )
    call output ( '==========================================', advance='yes' )
    call output ( 'subsection name ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'percent of total retrieval time', advance='yes' )
    total = sum(section_timings(1+num_section_times:))
    do elem = 1, num_retrieval_times
      call GetStringElement(retrieval_names, section_name, elem, countEmpty)
      percent = 100 * section_timings(num_section_times+elem) / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( section_timings(num_section_times+elem), FORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( percent, FORMAT=PCTFORM, advance='yes' )
    enddo
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 7, advance='no' )
    call output ( final-total, FORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( 100*(final-total)/final, FORMAT=PCTFORM, advance='yes' )
  end subroutine dump_section_timings

!=============================================================================
END MODULE MLSL2Timings
!=============================================================================

!
! $Log$
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
