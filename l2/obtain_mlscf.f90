! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module OBTAIN_MLSCF

! Open and close the MLSCF
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use SDPToolkit, only: Pgs_io_gen_closeF, Pgs_io_gen_openF, PGS_S_SUCCESS, &
    & PGSd_IO_Gen_RSeqFrm

  implicit NONE

  private
  public :: Close_MLSCF, Open_MLSCF

  ! =====  Private declarations  =======================================

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
     "$id: obtain_mlscf.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  CLOSE_MLSCF  -----
  subroutine CLOSE_MLSCF ( CF_Unit )

    integer, intent(in) :: CF_Unit

    character (LEN=32) :: MNEMONIC
    character (LEN=256) :: MSG
    integer :: RETURN_STATUS

    return_Status = Pgs_io_gen_closeF ( CF_Unit )

    if ( return_Status /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( return_Status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Error closing L2CF:  '//mnemonic//' '//msg)
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

    L2CF_Version = 1
    return_Status = Pgs_io_gen_openF ( MLSPCF_Start, PGSd_IO_Gen_RSeqFrm, &
                                      0, CF_Unit, L2CF_Version)

    if ( return_Status /= PGS_S_SUCCESS ) then

      call Pgs_smf_getMsg ( return_Status, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Error opening MLSCF:  "//mnemonic//"  "//msg)
    end if

  end subroutine OPEN_MLSCF

end module OBTAIN_MLSCF

! $Log$
! Revision 2.2  2001/02/28 01:58:07  vsnyder
! Get pcf # and unit # from arguments
!
