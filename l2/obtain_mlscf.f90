! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module OBTAIN_MLSCF

! Open and close the MLSCF
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSFiles, only: MLS_io_gen_openf, MLS_io_gen_closef
  USE output_m, only: output
  use SDPToolkit, only: PGS_S_SUCCESS, PGSd_IO_Gen_RSeqFrm
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF

  implicit NONE

  private
  public :: Close_MLSCF, Open_MLSCF

  ! =====  Private declarations  =======================================

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
     "$id: obtain_mlscf.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  integer, private :: ERROR

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  CLOSE_MLSCF  -----
  subroutine CLOSE_MLSCF ( CF_Unit, return_status )

    integer, intent(in) :: CF_Unit
    integer, intent(out) :: return_status

    character (LEN=32) :: MNEMONIC
    character (LEN=256) :: MSG
!    integer :: RETURN_STATUS

    error = 0
!    return_Status = Pgs_io_gen_closeF ( CF_Unit )

    return_Status = Mls_io_gen_closeF ( 'pg', CF_Unit )

    if ( return_Status /= PGS_S_SUCCESS ) then
		call announce_error(0, 'Error closing L2CF', &
      & error_number=return_Status)
    end if

  end subroutine CLOSE_MLSCF

  ! -------------------------------------------------  OPEN_MLSCF  -----
  subroutine OPEN_MLSCF ( MLSPCF_Start, CF_Unit, return_status )

    integer, intent(in) :: MLSPCF_Start
    integer, intent(out) :: CF_Unit
    integer, intent(out) :: return_status

    integer :: L2CF_VERSION
    integer :: record_length
    character (LEN=32) :: MNEMONIC
    character (LEN=256) :: MSG

    error = 0
    L2CF_Version = 1
!    return_Status = Pgs_io_gen_openF ( MLSPCF_Start, PGSd_IO_Gen_RSeqFrm, &
!                                      0, CF_Unit, L2CF_Version)

    CF_Unit = Mls_io_gen_openF ( 'pg', .true., return_Status, record_length, &
      & PGSd_IO_Gen_RSeqFrm, &
      & thePC=MLSPCF_Start)

    if ( return_Status /= PGS_S_SUCCESS ) then
		call announce_error(0, "Error opening MLSCF", &
      & error_number=return_Status)
    end if

  end subroutine OPEN_MLSCF

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
  & error_number )
  
   ! Arguments
	
    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number

    ! Local
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.
 
      error = 0
    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
 
    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error: ", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
        call output ( 'Error number ' )
        call output ( error_number, advance='yes' )
      end if
    end if

!===========================
  end subroutine announce_error
!===========================

end module OBTAIN_MLSCF

! $Log$
! Revision 2.8  2001/07/18 23:57:57  pwagner
! Returns error from close_mlscf
!
! Revision 2.7  2001/05/02 23:33:48  pwagner
! Replace pgs_io_gen.. routines with mls_..
!
! Revision 2.6  2001/04/16 23:43:17  pwagner
! Returns returnStatus
!
! Revision 2.5  2001/04/12 22:19:33  vsnyder
! Improved an error message
!
! Revision 2.4  2001/04/05 23:42:10  pwagner
! Added announce_error, deleted all MLSMessages
!
! Revision 2.3  2001/02/28 02:02:21  vsnyder
! Remove unused reference to MLSPCF2
!
! Revision 2.2  2001/02/28 01:58:07  vsnyder
! Get pcf # and unit # from arguments
!
