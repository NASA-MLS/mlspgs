! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!===================================================================
module ObtainClimatology 

! Provides subroutines to access climatology files in L3 ascii format
!===================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dates_module, only: CCSDS2TAI
  use GriddedData, only: DestroyGridTemplateContents, GriddedData_T, &
    & AddGridTemplateToDatabase, l3ascii_read_field, &
    & make_log_axis
  use MLSCommon, only: LineLen, R8
  use MLSMessagemodule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSPCF2, only: mlspcf_l2clim_end, mlspcf_l2clim_start
  use MLSStrings, only: Capitalize, Count_Words, ReadCompleteLineWithoutComments
  use SDPToolkit, only: Pgs_io_gen_openF, PGS_S_SUCCESS, PGSd_IO_Gen_RSeqFrm
! use VerticalCoordinate

  implicit none
  private
  public :: OBTAIN_CLIM, READ_CLIMATOLOGY

  private :: Id, moduleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = & 
     "$Id$"
  character(len=*), parameter :: moduleName="$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  OBTAIN_CLIM  -----
  !=====================================================================
  subroutine OBTAIN_CLIM ( aprioriData, root )
  !=====================================================================
    !Arguments 
    type (GriddedData_T), dimension(:), pointer :: aprioriData 
    ! Input a priori database
    integer, intent(in) :: root        ! Root of the L2CF abstract syntax tree

    !Local Variables

    type (GriddedData_T):: qty
    character (LEN=256) :: msg, mnemonic
    integer:: CliUnit, processCli, returnStatus, version

    logical :: end_of_file = .FALSE.

    do CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

!     Open one Climatology file as a generic file for reading
      version = 1
      returnStatus = Pgs_io_gen_openF ( CliUnit, PGSd_IO_Gen_RSeqFrm, 0, &
                                        processCli, version )
      if ( returnStatus /= PGS_S_SUCCESS ) then

        call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Error opening Climatology file:  "//mnemonic//" "//msg)

      end if


      do while (.NOT. end_of_file)

        call l3ascii_read_field ( processCli, qty, end_of_file)
        returnStatus = AddGridTemplateToDatabase(aprioriData, qty)
        call DestroyGridTemplateContents ( qty )

      end do !(.not. end_of_file)

    end do ! CliUnit = mlspcf_l2clim_start, mlspcf_l2clim_end

  return
  !============================
  end subroutine OBTAIN_CLIM
  !============================

  ! --------------------------------------------------  READ_NCEP  -----
  SUBROUTINE READ_CLIMATOLOGY ( fname, data_array )
  ! --------------------------------------------------
  ! Brief description of program
  ! This subroutine reads a NCEP correlative file and returns
  ! the data_array to the caller

  ! Arguments

  character*(*), intent(in) :: fname			! Physical file name
  real(R8) ::  data_array(:,:,:)

	END SUBROUTINE READ_CLIMATOLOGY

! =====     Private Procedures     =====================================

!============================
end module ObtainClimatology
!============================
! $Log$
! Revision 2.3  2001/03/03 00:11:29  pwagner
! Began transformations to act like L2GPData module for Gridded data
!
! Revision 2.2  2001/02/21 00:37:51  pwagner
! Uses more of GriddedData
!
! Revision 2.1  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.0  2000/09/05 18:57:06  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

