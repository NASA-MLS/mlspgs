!=================================
PROGRAM PCNumberTest ! tests MLSFiles module
!=================================

   USE MLSFiles , ONLY: GetPCFromRef
   USE MLSStrings , ONLY: Capitalize

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the GetPCFromRef subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, enter a molecule name; e.g., "H2O" (w/o the "'s)

! Variables

	LOGICAL, PARAMETER                :: caseSensitive = .false.
   character(LEN=32)                 :: theMolecule
   character(LEN=132)                :: FileName
   integer                           :: returnStatus, ErrType

	DO

! Prompt for input

 	  print *, 'Enter a molecule name'
		read(*, '(A32)') theMolecule
      print *, 'You entered ', trim(theMolecule)
		IF(theMolecule == ' ' .or. Capitalize(theMolecule(1:4)) == 'STOP') THEN
   		print *, 'GAME OVER'
			STOP
		ENDIF
      
!	Process theMolecule
	   returnStatus = GetPCFromRef(theMolecule, 60000, 60018, &
       & caseSensitive, ErrType, ExactName=FileName)
      print *, 'returnStatus: ', returnStatus
      print *, 'ErrType: ', ErrType
      print *, 'ExactName: ', trim(FileName)
	ENDDO
!==================
END PROGRAM PCNumberTest
!==================
integer function Pgs_fake_getReference(arg, version, fileName)
! arguments
  integer, intent(in)                  :: arg
  integer, intent(in)                  :: version
  character(LEN=*), intent(out)        :: fileName
  ! Internal variables
  integer, parameter                   :: PCOFFSET = 60000
  integer, parameter                   :: FullNameLen = 132
  character(LEN=*), dimension(17), parameter :: ALL_FILES = (/ &    
    & '/my/path/to/here/MLS_AURA_L2GP_TEMP_V1_DATE_WHEN_TAKEN.dat', &       
    & '/my/path/to/here/MLS_AURA_L2GP_GPH_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_H2O_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_HNO3_V1_DATE_WHEN_TAKEN.dat', &       
    & '/my/path/to/here/MLS_AURA_L2GP_O3_V1_DATE_WHEN_TAKEN.dat  ', &       
    & '/my/path/to/here/MLS_AURA_L2GP_HCL_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_CLO_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_N2O_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_CO_V1_DATE_WHEN_TAKEN.dat  ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_OH_V1_DATE_WHEN_TAKEN.dat  ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_RHI_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_SO2_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_HO2_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_BRO_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_HOCL_V1_DATE_WHEN_TAKEN.dat', &      
    & '/my/path/to/here/MLS_AURA_L2GP_HCN_V1_DATE_WHEN_TAKEN.dat ', &      
    & '/my/path/to/here/MLS_AURA_L2GP_BULL_V1_DATE_WHEN_TAKEN.dat'/) 
    
  ! Executable statements
  select case (arg-PCOFFSET)
    case(:0)
      Pgs_fake_getReference = -1
    case(1:size(ALL_FILES))
      Pgs_fake_getReference = 0
      fileName = ALL_FILES(arg-PCOFFSET)
    case default
      Pgs_fake_getReference = 1
  end select

     
end function Pgs_fake_getReference
!# $Log$
