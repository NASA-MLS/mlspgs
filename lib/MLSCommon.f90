! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  USE SDPTOOLKIT, only: PGSd_PC_FILE_PATH_MAX

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module contains simple definitions which are common to all the MLS PGS
  ! f90 software.

  ! Firstly, these are standard numerical types, copied from HCP
  ! (again with my change in case, sorry Hugh!)

  INTEGER, PARAMETER:: i1=SELECTED_INT_KIND(2)
  INTEGER, PARAMETER:: i2=SELECTED_INT_KIND(4)
  INTEGER, PARAMETER:: i4=SELECTED_INT_KIND(7)
  INTEGER, PARAMETER:: r4=SELECTED_REAL_KIND(5)
  INTEGER, PARAMETER:: r8=SELECTED_REAL_KIND(13)

  ! Now choose the precision we want by preference (may automate this through
  ! make later on, with perl or m4 or something).
  INTEGER, PARAMETER:: rp=r8
  INTEGER, PARAMETER:: ip=i4

  ! Now we have the lengths for various strings

  INTEGER, PARAMETER :: NameLen=32
  INTEGER, PARAMETER :: LineLen=132

! Shouldn't this be PGSd_PC_FILE_PATH_MAX ?
  INTEGER, PARAMETER :: FileNameLen=max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132

  ! --------------------------------------------------------------------------
  
  ! The next datatype describes the information on the L1B data files in use

  TYPE L1BInfo_T
    INTEGER :: L1BOAId=0     ! The HDF ID (handle) for the L1BOA file
    ! Id(s) for the L1BRAD file(s)
    INTEGER, DIMENSION(:), POINTER :: L1BRADIds=>NULL()
    CHARACTER (LEN=FileNameLen) :: L1BOAFileName=""  ! L1BOA file name
    CHARACTER (LEN=FileNameLen), DIMENSION(:), POINTER :: &
         & L1BRADFileNames=>NULL()
  END TYPE L1BInfo_T

  ! --------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  TYPE MLSChunk_T
    INTEGER :: firstMAFIndex   ! Index of first MAF in the chunk
    INTEGER :: lastMAFIndex    ! Index of last MAF in the chunk
    INTEGER :: noMAFsLowerOverlap ! Number of MAFs in the lower overlap region
    INTEGER :: noMAFsUpperOverlap ! Number of MAFs in the upper overlap region
    INTEGER :: accumulatedMAFs ! Number of non overlapped MAFs before this.
  END TYPE MLSChunk_T

  ! --------------------------------------------------------------------------

  ! The TAI93 time range

  TYPE TAI93_Range_T
    REAL(r8) :: startTime ! TAI93 format
    REAL(r8) :: endTime   ! TAI93 format
  END TYPE TAI93_Range_T
  ! --------------------------------------------------------------------------

  contains

  ! -------------------------------------------- FindFirst --------------
  integer function FindFirst ( condition )
    ! Find the first logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION

    ! Local variables
    integer :: I                        ! Loop counter

    ! Executable code
    FindFirst = 0
    do i = 1, size(condition)
      if ( condition(i) ) then
        FindFirst = i
        return
      end if
    end do
  end function FindFirst


!=============================================================================
END MODULE MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.9  2002/01/09 23:51:27  pwagner
! Connected FileNameLen with PGSd_PC_FILE_PATH_MAX
!
! Revision 2.8  2001/11/14 18:03:32  livesey
! Changed FindFirst to return 0 not -1 if not found
!
! Revision 2.7  2001/09/09 02:47:58  livesey
! Moved FindFirst into MLSCommon
!
! Revision 2.6.2.3  2001/09/09 01:53:27  livesey
! Bug fix
!
! Revision 2.6.2.2  2001/09/09 01:35:46  livesey
! Moved FindFirst in from MLSL2Common
!
! Revision 2.6.2.1  2001/09/08 22:32:24  livesey
! Added RP and IP
!
! Revision 2.6  2001/04/20 23:10:53  livesey
! Initialised parameters in L1BINFO
!
! Revision 2.5  2001/03/10 18:48:17  livesey
! Really nullified the pointer!
!
! Revision 2.4  2001/03/10 07:06:46  livesey
! Nullified L1BRadfileNames in L1BInfo
!
! Revision 2.3  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.2  2001/01/26 23:46:35  pwagner
! Restored L1BInfo from l1/MLSL1Common back to lib/MLSCommon
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
