! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
!MODULE SDPToolkitSubstitute     ! Substitute for the essential toolkit routines
!=============================================================================

!  IMPLICIT NONE

!  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
!  CHARACTER (LEN=256) :: Id = &
!       "$Id$"
!  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module gives substitutes for the essential toolkit routines used by
  ! library code used in both the toolkit and non toolkit environment.

  ! Some of the functions don't work. Code that calls them will need
  ! to be linked to the _real_ SDP Toolkit library PGSTK.a. They are in here
  ! to satisfy some F90 compilers which check that a routine exists if it
  ! has an interface block defined in another module, even if it is not
  ! called in a particular case of USEing that module. There shoulld be a 
  ! routine in here for every SDP routine that has an interface in 
  ! SDPToolkit.f90

!CONTAINS

       INTEGER FUNCTION PGS_SMF_GenerateStatusReport(msg)
             CHARACTER (LEN=*):: msg
             !INTEGER :: PGS_SMF_GenerateStatusReport

             PRINT*,msg
             PGS_SMF_GenerateStatusReport=0
       END FUNCTION PGS_SMF_GenerateStatusReport

      INTEGER FUNCTION PGS_IO_Gen_OpenF(file_logical,file_access, &
           & record_length, file_handle, file_version)
        INTEGER, INTENT(IN) :: file_logical
        INTEGER, INTENT(IN) :: file_access
        INTEGER, INTENT(IN) :: record_length
        INTEGER, INTENT(OUT) :: file_handle
        INTEGER, INTENT(IN) :: file_version
        
        PGS_IO_Gen_OpenF=-99
      END FUNCTION PGS_IO_Gen_OpenF


      INTEGER FUNCTION PGS_IO_Gen_CloseF(file_handle)
        INTEGER, INTENT(IN) :: file_handle
        PGS_IO_Gen_CloseF=-99
      END FUNCTION PGS_IO_Gen_CloseF

      INTEGER FUNCTION Pgs_pc_getReference(file_handle, &
         file_version, physicalfilename)
        INTEGER, INTENT(IN) :: file_handle
        INTEGER, INTENT(INOUT) :: file_version
        character (LEN=*), INTENT(OUT) :: physicalfilename
        Pgs_pc_getReference=-99
      END FUNCTION Pgs_pc_getReference


      subroutine Pgs_smf_getMsg ( CODE, MNEMONIC, MSG )
        integer, intent(out) :: CODE              ! Previously stored code
        character(len=*), intent(out) :: MNEMONIC ! Previously stored mnemonic
        character(len=*), intent(out) :: MSG      ! Previously stored message
      end subroutine Pgs_smf_getMsg

      INTEGER FUNCTION PGS_TD_TAItoUTC(sectai93,asciiutc)
        DOUBLE PRECISION, INTENT(IN) :: sectai93
        CHARACTER(LEN=27), INTENT(OUT) :: asciiutc
        asciiutc="Use the Real PGSTK"
        PGS_TD_TAItoUTC=-99
      END FUNCTION PGS_TD_TAItoUTC

!=============================================================================
!END MODULE SDPToolkitSubstitute
!=============================================================================

!
! $Log$
! Revision 2.2  2000/11/08 10:21:04  pumphrey
! Added return values to functions to prevent compile errors
!
! Revision 2.1  2000/11/06 18:52:22  pumphrey
! Added routines so that every routine in SDPToolkit.f90 is represented.
! The routines do nothing, they are just to satisfy the link stage with
! some of the odder Fortran 9x compilers (NA FortranPlus in fact)
!
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.3  1999/12/03 00:14:18  livesey
! Nightly checkin
!
