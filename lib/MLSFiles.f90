! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE MLSFiles               ! Utility file routines
!===============================================================================
   USE MLSCommon, only: i4
   USE MLSStrings, only: Capitalize
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

  CONTAINS

  ! ---------------------------------------------  GetPCFromRef  -----

  ! This function takes a FileName as an arg and a range of PC numbers
  ! [PCBottom, PCTop] which are integers
  ! It returns thePC corresponding to the FileName
  ! If the FileName is not found, it sets ErrType=NAMENOTFOUND
  ! otherwise ErrType=0
  
  ! This is useful because all the Toolbox routines refer to files
  ! by their PC numbers, not their names
  
  FUNCTION GetPCFromRef(FileName, PCBottom, PCTop, &
  & caseSensitive, ErrType, versionNum) RESULT (thePC)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN)   :: FileName
    INTEGER(i4),  INTENT(IN)   :: PCBottom, PCTop
    INTEGER(i4)                :: thePC
    INTEGER(i4),  INTENT(OUT)  :: ErrType
    LOGICAL,  INTENT(IN)       :: caseSensitive
    INTEGER(i4),  OPTIONAL     :: versionNum
  
    ! Local variables
	INTEGER, PARAMETER :: NAMENOTFOUND=-1
	INTEGER, PARAMETER :: INVALIDPCRANGE=NAMENOTFOUND-1
	INTEGER, PARAMETER :: MAXNAMELENGTH=132
	
	CHARACTER (LEN=MAXNAMELENGTH) :: MatchName, TryName
	INTEGER                       :: version, returnStatus
		
	!
	thePC = 0
	IF(PCTop < PCBottom) THEN
		ErrType = INVALIDPCRANGE
		RETURN
	ENDIF
	
	IF(caseSensitive) THEN
		MatchName = FileName
	ELSE
		MatchName = Capitalize(FileName)
	ENDIF

	IF(PRESENT(versionNum)) THEN
		version = versionNum
	ELSE
		version = 1
	ENDIF

	ErrType = NAMENOTFOUND

	DO thePC = PCBottom, PCTop

            returnStatus = Pgs_pc_getReference(thePC, version, &
              & TryName)
            
            if ( returnStatus == PGS_S_SUCCESS ) then

		IF(.NOT. caseSensitive) THEN
			TryName = Capitalize(TryName)
		ENDIF

              if ( INDEX(TryName, TRIM(MatchName)) /= 0 )then
					ErrType = 0
                EXIT
					endif

				endif

	ENDDO

  END FUNCTION GetPCFromRef
	
!====================
END MODULE MLSFiles
!====================

!
! $Log$
! Revision 2.1  2001/03/07 01:02:37  pwagner
! First commit
!
!
