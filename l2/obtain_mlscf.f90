! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module OBTAIN_MLSCF

! Open and close the MLSCF
  use LEXER_CORE, only: PRINT_SOURCE
!  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  USE output_m, only: output
  use SDPToolkit, only: Pgs_io_gen_closeF, Pgs_io_gen_openF, PGS_S_SUCCESS, &
    & PGSd_IO_Gen_RSeqFrm
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
  subroutine CLOSE_MLSCF ( CF_Unit )

    integer, intent(in) :: CF_Unit

    character (LEN=32) :: MNEMONIC
    character (LEN=256) :: MSG
    integer :: RETURN_STATUS

    error = 0
    return_Status = Pgs_io_gen_closeF ( CF_Unit )

    if ( return_Status /= PGS_S_SUCCESS ) then
!      call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!      call MLSMessage ( MLSMSG_Error, ModuleName, &
!        & 'Error closing L2CF:  '//mnemonic//' '//msg)
		call announce_error(0, 'Error closing L2CF')
    end if

  end subroutine CLOSE_MLSCF

  ! -------------------------------------------------  OPEN_MLSCF  -----
  subroutine OPEN_MLSCF ( MLSPCF_Start, CF_Unit )

    integer, intent(in) :: MLSPCF_Start
    integer, intent(out) :: CF_Unit

    integer :: L2CF_VERSION
    character (LEN=32) :: MNEMONIC
    character (LEN=256) :: MSG
    integer :: RETURN_STATUS

    error = 0
    L2CF_Version = 1
    return_Status = Pgs_io_gen_openF ( MLSPCF_Start, PGSd_IO_Gen_RSeqFrm, &
                                      0, CF_Unit, L2CF_Version)

    if ( return_Status /= PGS_S_SUCCESS ) then
!      call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!      call MLSMessage ( MLSMSG_Error, ModuleName, &
!        & "Error opening MLSCF:  "//mnemonic//"  "//msg)
		call announce_error(0, "Error opening MLSCF")
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
	if(present(use_toolkit)) then
		just_print_it = use_toolkit
	elseif(default_output_by_toolkit) then
		just_print_it = .false.
	else
		just_print_it = .true.
	endif
	
	if(.not. just_print_it) then
    error = max(error,1)
    call output ( '***** At ' )

	if(lcf_where > 0) then
	    call print_source ( source_ref(lcf_where) )
		else
    call output ( '(no lcf node available)' )
		endif

    call output ( ': ' )
    call output ( "The " );
	if(lcf_where > 0) then
    call dump_tree_node ( lcf_where, 0 )
		else
    call output ( '(no lcf tree available)' )
		endif

		CALL output("Caused the following error:", advance='yes', &
		& from_where=ModuleName)
		CALL output(trim(full_message), advance='yes', &
		& from_where=ModuleName)
		if(present(error_number)) then
			CALL output('error number ', advance='no')
			CALL output(error_number, places=9, advance='yes')
		endif
	else
		print*, '***Error in module ', ModuleName
		print*, trim(full_message)
		if(present(error_number)) then
			print*, 'error number ', error_number
		endif
	endif

!===========================
  end subroutine announce_error
!===========================

end module OBTAIN_MLSCF

! $Log$
! Revision 2.3  2001/02/28 02:02:21  vsnyder
! Remove unused reference to MLSPCF2
!
! Revision 2.2  2001/02/28 01:58:07  vsnyder
! Get pcf # and unit # from arguments
!
