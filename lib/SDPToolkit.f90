! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE SDPToolkit               ! F90 interface to SDP Toolkit.
!===============================================================================
   USE MLSCommon
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

! Contents:

! PGS_CSC_4.f
! PGS_EPH_5.f
! PGS_IO.f
! PGS_IO_1.f
! PGS_MEM_7.f
! PGS_PC.f
! PGS_PC_9.f
! PGS_SMF.f
! PGS_TD.f
! PGS_TD_3.f

! Plust interfaces for toolkit routines.

! Remarks:  This module contains include files required to use the toolkit.
! Note that the include files with numbers in the names are generated 
! automagically when you install the toolkit. They are *not* in the toolkit
! tarball. If the toolkit has not installed properly on your system, those
! files will not be present and hence you will find that this file will 
! not compile

   INCLUDE 'PGS_CBP.f'
   INCLUDE 'PGS_CBP_6.f'
   INCLUDE 'PGS_CSC_4.f'
   INCLUDE 'PGS_EPH_5.f'
   INCLUDE 'PGS_IO.f'
   INCLUDE 'PGS_IO_1.f'
   INCLUDE 'PGS_MEM_7.f'
   INCLUDE 'PGS_PC.f'
   INCLUDE 'PGS_PC_9.f'
   INCLUDE 'PGS_SMF.f'
   INCLUDE 'PGS_TD.f'
   INCLUDE 'PGS_TD_3.f'
   INCLUDE 'PGS_MET_13.f'
   INCLUDE 'PGS_MET.f'


! Now define f90 interfaces for some toolkit routines.

   INTERFACE
      INTEGER FUNCTION PGS_IO_Gen_OpenF(file_logical,file_access, &
           & record_length, file_handle, file_version)
        INTEGER, INTENT(IN) :: file_logical
        INTEGER, INTENT(IN) :: file_access
        INTEGER, INTENT(IN) :: record_length
        INTEGER, INTENT(OUT) :: file_handle
        INTEGER, INTENT(IN) :: file_version
      END FUNCTION PGS_IO_Gen_OpenF

      INTEGER FUNCTION PGS_IO_Gen_CloseF(file_handle)
        INTEGER, INTENT(IN) :: file_handle
      END FUNCTION PGS_IO_Gen_CloseF

      INTEGER FUNCTION Pgs_pc_getReference(file_handle, file_version, physicalfilename)
        INTEGER, INTENT(IN) :: file_handle
        INTEGER, INTENT(INOUT) :: file_version
        character (LEN=*), INTENT(OUT) :: physicalfilename
      END FUNCTION Pgs_pc_getReference

      INTEGER FUNCTION PGS_SMF_GenerateStatusReport(msg)
        CHARACTER (LEN=*) :: msg
      END FUNCTION PGS_SMF_GenerateStatusReport

      subroutine Pgs_smf_getMsg ( CODE, MNEMONIC, MSG )
        integer, intent(out) :: CODE              ! Previously stored code
        character(len=*), intent(out) :: MNEMONIC ! Previously stored mnemonic
        character(len=*), intent(out) :: MSG      ! Previously stored message
      end subroutine Pgs_smf_getMsg

      INTEGER FUNCTION PGS_TD_TAItoUTC(sectai93,asciiutc)
        DOUBLE PRECISION, INTENT(IN) :: sectai93
        CHARACTER(LEN=27), INTENT(OUT) :: asciiutc
      END FUNCTION PGS_TD_TAItoUTC

   END INTERFACE

!   INTEGER, EXTERNAL :: PGS_SMF_GenerateStatusReport
!   INTEGER, EXTERNAL :: PGS_IO_Gen_OpenF

! Parameters used by multiple modules

!  CHARACTER (LEN=*), PARAMETER :: earthModel = 'WGS84'
   CHARACTER (LEN=*), PARAMETER :: earthModel = 'GRS-80'

   INTEGER, PARAMETER :: max_orbits = 16
   INTEGER, PARAMETER :: spacecraftId = 6666

   REAL, PARAMETER :: Deg2Rad = 1.7453293E-02
   REAL, PARAMETER :: PI = 3.141592654
   REAL, PARAMETER :: Rad2Deg = 57.295780

!====================
END MODULE SDPToolkit
!====================

!
! $Log$
! Revision 2.3  2000/09/27 15:01:08  pumphrey
! Reinstated "numbered" Toolkit include files now I understand
! where they come from. (Duh. Groan.)
!
! Revision 2.2  2000/09/26 14:18:02  pumphrey
! Changed an arg of PGS_PC_Getreference to INOUT.
!
! Revision 2.1  2000/09/19 11:16:38  pumphrey
! Removed INCLUDEs of files that are not supplied with the SDP toolkit
!
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 02:00:46  vsnyder
! Cosmetic changes
!
