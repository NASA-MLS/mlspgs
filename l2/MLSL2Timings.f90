! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
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
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! This module simply contains initial settings and accumulated values.

  ! The following public settings are stored here; they may be set by MLSL2
  logical          :: SECTION_TIMES = .false.  ! Show times in each section
  logical          :: TOTAL_TIMES = .false.    ! Show total times from start

  logical, private   :: COUNTEMPTY = .false.     ! Any sections named ' '?
  logical, private   :: ALLOWUNKNOWNNAMES = .false.  ! Any unknown section names

  character*(*), parameter           :: section_names = &
    & 'main,open_init,global_settings,signals,spectroscopy,' // &
    & 'read_apriori,scan_divide,construct,fill,retrieve,join,output'
  integer, parameter                 :: num_section_times = 12
  character*(*), parameter           :: retrieval_names = &
    & 'newton_solver,cholesky_factor,cholesky_solver,cholesky_invert,' // &
    & 'full_fwm,fullcloud_fwm,scan_fwm,twod_scan_fwm,linear_fwm,' // &
    & 'low_cloud,high_cloud,sids,form_normeq,tikh_reg'
  integer, parameter                 :: num_retrieval_times = 14
  ! dimension of the following is +2 to allow possible unknown
  ! section names and unknown retrieval names
  real, dimension(num_section_times+num_retrieval_times+2), &
    & save                           :: section_timings = 0.
  integer, parameter                 :: unknown_section = &
    &                               num_section_times+num_retrieval_times + 1
  integer, parameter                 :: unknown_retrieval = unknown_section + 1

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------  add_to_retrieval_timing  -----
  subroutine add_to_retrieval_timing( section_name, t1 )
  ! Add current elapsed retrieval section time to total so far for section_name

  ! Formal arguments
    character(LEN=*), intent(in):: section_name   ! One of the retr. sect_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    real                        :: t2
    real, save                  :: myLastTime

  ! Executable
      if ( present(t1) ) myLastTime = t1

      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
        if ( ALLOWUNKNOWNNAMES ) then
          call time_now ( t2 )
          section_timings(unknown_retrieval) = &
            & section_timings(unknown_retrieval) + t2 - myLastTime
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find retrieval section name ' // section_name // &
        & ' among list ' // retrieval_names )
        endif
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
    real, save                  :: myLastTime

  ! Executable
    if ( present(t1) ) myLastTime = t1
    elem = StringElementNum(section_names, LowerCase(section_name), countEmpty)
    if ( elem < 1 .or. elem > num_section_times ) then
      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
        if ( ALLOWUNKNOWNNAMES ) then
          call time_now ( t2 )
          section_timings(unknown_section) = &
            & section_timings(unknown_section) + t2 - myLastTime
        else
          call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Unable to find section name ' // section_name // &
          & ' among list ' // section_names // ',' // retrieval_names )
        endif
      else
        call time_now ( t2 )
        section_timings(num_section_times+elem) = &
          & section_timings(num_section_times+elem) + t2 - myLastTime
      endif
    else
      call time_now ( t2 )
      section_timings(elem) = section_timings(elem) + t2 - myLastTime
    endif
    if ( present(t1) ) call time_now ( t1 )
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
    real                            :: elem_time
    character(LEN=LEN(TIMEFORMBIG)) :: TIMEFORM
    logical                         :: Unknown_nonzero
    integer                         :: num_elems

  ! Executable
    Unknown_nonzero = (section_timings(unknown_section) > 0.)
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Level 2 section timings : ', advance='yes' )
    call output ( '==========================================', advance='yes' )
    call output ( 'section name ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'percent ', advance='yes' )
    total = sum(section_timings(1:num_section_times)) + &
     & section_timings(unknown_section)
    call time_now ( final )
    final = max(final, total)
    if ( final <= 0.0 ) final = 1.0       ! Just so we don't divide by 0
    if ( final < 0.5 .or. final > 99999.99 ) then
      TIMEFORM = TIMEFORMBIG
    else
      TIMEFORM = TIMEFORMSMALL
    endif
    if ( Unknown_nonzero ) then
      num_elems = num_section_times + 1
    else
      num_elems = num_section_times
    endif
    do elem = 1, num_section_times
      if ( elem > num_section_times ) then
        elem_time = section_timings(unknown_section)
        section_name = 'unknown name'
      else
        elem_time = section_timings(elem)
        call GetStringElement(section_names, section_name, elem, countEmpty)
      endif
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time, FORMAT=TIMEFORM, advance='no' )
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
    total = sum(section_timings(1+num_section_times:)) - &
      & section_timings(unknown_section)
    Unknown_nonzero = (section_timings(unknown_retrieval) > 0.)
    if ( Unknown_nonzero ) then
      num_elems = num_retrieval_times + 1
    else
      num_elems = num_retrieval_times
    endif
    do elem = 1, num_elems
      if ( elem > num_retrieval_times ) then
        elem_time = section_timings(unknown_retrieval)
        section_name = 'unknown name'
      else
        elem_time = section_timings(num_section_times+elem)
        call GetStringElement(retrieval_names, section_name, elem, countEmpty)
      endif
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time, FORMAT=TIMEFORM, advance='no' )
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
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END MODULE MLSL2Timings
!=============================================================================

!
! $Log$
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
