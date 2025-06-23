! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TABLE_DUMPER
  use DECLARATION_TABLE, only: ENUM_VALUE, EXPRN, LABEL, LOG_VALUE, &
                               NUM_VALUE, RANGE, STR_RANGE, STR_VALUE, &
                               TYPE_NAMES
  use INTRINSIC, only: PHYQ_INDICES
  use MLSCF, only: Mlscf_T, MlscfCell_T
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING
  implicit NONE

  private
  public :: DUMP_TABLE

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
  subroutine DUMP_TABLE ( L2CF_DATA )
    type(mlscf_t), intent(in) :: L2CF_DATA   ! The L2CF data (duh)
    integer :: I, J, K                       ! Loop inductors
    call output ( l2cf_data%NoSections )
    call output ( ' Sections', advance='yes' )
    do i = 1, l2cf_data%NoSections
      call output ( 'Section: ' )
      call output ( trim(l2cf_data%sections(i)%MlscfSectionName), &
                    advance='yes' )
      if ( l2cf_data%sections(i)%NoSectionDefs > 0 ) then
        call output ( ' Definitions (' )
        call output (l2cf_data%sections(i)%NoSectionDefs )
        call output ( '):', advance='yes' )
      end if
      do j = 1, l2cf_data%sections(i)%NoSectionDefs
        call dump_cell ( l2cf_data%sections(i)%Cells(j), '  ' )
      end do
      if ( l2cf_data%sections(i)%NoSectionEntries > 0 ) then
        call output ( ' Specifications (' )
        call output ( l2cf_data%sections(i)%NoSectionEntries )
        call output ( '):', advance='yes' )
      end if
      do j = 1, l2cf_data%sections(i)%NoSectionEntries
        call output ( '  Specification: ' )
        if ( l2cf_data%sections(i)%Entries(j)%MlscfLabelName /= ' ' ) then
          call output ( trim(l2cf_data%sections(i)%Entries(j)%MlscfLabelName) )
          call output ( ': ' )
        end if
        call output ( trim(l2cf_data%sections(i)%Entries(j)%MlscfEntryName), &
                      advance = 'yes' )
        do k = 1, l2cf_data%sections(i)%Entries(j)%MlscfEntryNoKeys
          call dump_cell ( l2cf_data%sections(i)%Entries(j)%cells(k), '   ' )
        end do
      end do
    end do
  end subroutine DUMP_TABLE

  subroutine DUMP_CELL ( CELL, PREFIX )
    type(MlscfCell_T), intent(in) :: CELL
    character(len=*) :: PREFIX
    call output ( prefix )
    call output ( trim(cell%Keyword) )
    call output ( ' = ' )
    select case ( cell%type )
    case ( exprn )     ! unevaluated
      call output ( 'Unevaluated expression' )
    case ( label )     ! A label reference
      call output ( trim(cell%charValue) )
      if ( cell%charRangeUpperBound /= ' ' ) then
        call output ( '.' )
        call output ( trim(cell%charRangeUpperBound) )
      end if
    case ( log_value ) ! logical value
      call output ( cell%realValue )
    case ( num_value ) ! numerical value
      call output ( cell%realValue )
      call output ( ', units = ' )
      call display_string ( phyq_indices(cell%units) )
    case ( range )     ! numerical range
      call output ( cell%realValue )
      call output ( ' : ' )
      call output ( cell%RangeUpperBound )
      call output ( ', units = ' )
      call display_string ( phyq_indices(cell%units) )
    case ( str_range ) ! string range -- for dates
      call output ( trim(cell%charValue) )
      call output ( ' : ' )
      call output ( trim(cell%charRangeUpperBound) )
    case ( str_value, enum_value ) ! string value
      call output ( trim(cell%charValue) )
    case default
    end select
    call output ( ', TYPE: ' )
    call output ( trim(type_names(cell%type)) )
    if ( cell%more /= 0 ) then
      call output ( '; ' )
      call output ( cell%more )
      call output ( ' additional values' )
    end if
    call output ( '', advance='yes' )
  end subroutine DUMP_CELL
end module TABLE_DUMPER

! $Log$
! Revision 2.1  2000/10/11 18:57:28  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.2  2000/09/01 21:42:15  vsnyder
! Add "more" field to MLSCF_CELL type
!
