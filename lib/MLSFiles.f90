! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE MLSFiles               ! Utility file routines
!===============================================================================
   USE HDFEOS, only: gdopen, swopen
   USE machine, only: io_error
   USE MLSCommon, only: i4, NameLen
   USE MLSStrings, only: Capitalize, LowerCase
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg, &
  & PGSd_IO_Gen_RSeqFrm, PGSd_IO_Gen_RSeqUnf, & 
  & PGSd_IO_Gen_RDirFrm, PGSd_IO_Gen_RDirUnf, & 
  & PGSd_IO_Gen_WSeqFrm, PGSd_IO_Gen_WSeqUnf, & 
  & PGSd_IO_Gen_WDirFrm, PGSd_IO_Gen_WDirUnf, & 
  & PGSd_IO_Gen_USeqFrm, PGSd_IO_Gen_USeqUnf, & 
  & PGSd_IO_Gen_UDirFrm, PGSd_IO_Gen_UDirUnf, & 
  & PGSd_IO_Gen_ASeqFrm, PGSd_IO_Gen_ASeqUnf, PGS_IO_GEN_OpenF
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
!----------------------------------------------------------

  ! Now we have the legal unit numbers that files may be assigned

  INTEGER, PARAMETER :: bottom_unit_num=1
  INTEGER, PARAMETER :: top_unit_num=99

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
	INTEGER, PARAMETER :: MAXNAMELENGTH=NameLen
	
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
	
  ! ---------------------------------------------  mls_io_gen_openF  -----

  ! This function opens a generic file using either the toolbox
  ! or else a Fortran OPEN statement
  ! according to toolbox_mode
  ! It returns theFileHandle corresponding to the FileName or the PC

  ! If given a FileName as an arg and a range of PC numbers
  ! [PCBottom, PCTop] which are integers
  ! it will attempt to find a corresponding PC
  
  ! If given a PC it will attempt to find the corresponding FileName
  
  ! toolbox_mode                  meaning
  ! PGS_IO_GEN_OpenF              use PGS_IO_Gen_OpenF or fail
  !      swopen                   use swopen or fail
  !      gdopen                   use gdopen or fail
  !       open                    use Fortran or fail
  
  ! If the FileName is not found, it sets ErrType=NAMENOTFOUND
  ! otherwise ErrType=0
  
  ! This is useful because all the Toolbox routines refer to files
  ! by their PC numbers, not their names
  
  FUNCTION mls_io_gen_openF(toolbox_mode, caseSensitive, ErrType, &
  & record_length, FileAccessType, &
  & FileName, PCBottom, PCTop, versionNum, thePC) &
  &  RESULT (theFileHandle)

    ! Dummy arguments
    INTEGER(i4),  INTENT(OUT)  :: ErrType
    INTEGER(i4),  INTENT(OUT)  :: record_length
    LOGICAL,  INTENT(IN)       :: caseSensitive
    CHARACTER (LEN=*), INTENT(IN)   :: toolbox_mode
    INTEGER(i4), INTENT(IN)                :: FileAccessType
    INTEGER(i4)  :: theFileHandle
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN)   :: FileName
    INTEGER(i4),  OPTIONAL, INTENT(IN)   :: PCBottom, PCTop
    INTEGER(i4), OPTIONAL, INTENT(IN)                :: thePC
    INTEGER(i4),  OPTIONAL     :: versionNum
  
    ! Local variables

	INTEGER, PARAMETER :: UNKNOWNFILEACCESSTYPE=-999
	INTEGER, PARAMETER :: UNKNOWNTOOLBOXMODE=UNKNOWNFILEACCESSTYPE+1
	INTEGER, PARAMETER :: NOFREEUNITS=UNKNOWNTOOLBOXMODE+1
	INTEGER, PARAMETER :: MUSTSUPPLYFILENAMEORPC=NOFREEUNITS+1
	INTEGER, PARAMETER :: FH_ON_ERROR=-99
	INTEGER, PARAMETER :: KEYWORDLEN=12			! Length of keywords in OPEN(...)
	CHARACTER (LEN=NameLen) :: myName
	INTEGER(i4) :: myPC
	INTEGER                       :: version, returnStatus
   LOGICAL       :: tiedup
	CHARACTER (LEN=KEYWORDLEN) :: access, action, form, position, status
	INTEGER                       :: unit

	! begin

	! In case of premature return
	theFileHandle = FH_ON_ERROR

	IF(PRESENT(versionNum)) THEN
		version = versionNum
	ELSE
		version = 1
	ENDIF

	if(PRESENT(thePC)) then
		myPC = thePC
		returnStatus = Pgs_pc_getReference(thePC, version, &
              & myName)
	elseif(PRESENT(FileName)) then
		myName = FileName
		if(LowerCase(toolbox_mode(1:2)) == 'pg') then
			myPC = GetPCFromRef(FileName, PCBottom, PCTop, &
  &	 caseSensitive, returnStatus, versionNum)
  		endif
	else
		ErrType = MUSTSUPPLYFILENAMEORPC
		return
	ENDIF
	
	select case (LowerCase(toolbox_mode(1:2)))
	
	case('pg')
		if(returnStatus == 0) then
			ErrType = PGS_IO_Gen_OpenF(myPC, FileAccessType, record_length, &
			& theFileHandle, version)
		endif
	
	case('sw')
		if(returnStatus == 0) then
			theFileHandle = swopen(myName, FileAccessType)
			if(theFileHandle <= 0) then
				ErrType = min(theFileHandle, -1)
			else
				ErrType = 0
			endif
		else
			ErrType = returnStatus
		endif
	
	case('gd')
		if(returnStatus == 0) then
			theFileHandle = gdopen(myName, FileAccessType)
			if(theFileHandle <= 0) then
				ErrType = min(theFileHandle, -1)
			else
				ErrType = 0
			endif
		else
			ErrType = returnStatus
		endif
	
	case('op')
		theFileHandle = FH_ON_ERROR
		if(FileAccessType == PGSd_IO_Gen_RSeqFrm) then
			status = 'old'
			access = 'sequential'
			form = 'formatted'
			action = 'read'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_RSeqUnf) then
			status = 'old'
			access = 'sequential'
			form = 'unformatted'
			action = 'read'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_RDirFrm) then
			status = 'old'
			access = 'direct'
			form = 'formatted'
			action = 'read'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_RDirUnf) then
			status = 'old'
			access = 'direct'
			form = 'unformatted'
			action = 'read'
			position = 'rewind'

		elseif(FileAccessType == PGSd_IO_Gen_WSeqFrm) then
			status = 'new'
			access = 'sequential'
			form = 'formatted'
			action = 'write'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_WSeqUnf) then
			status = 'new'
			access = 'sequential'
			form = 'unformatted'
			action = 'write'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_WDirFrm) then
			status = 'new'
			access = 'direct'
			form = 'formatted'
			action = 'write'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_WDirUnf) then
			status = 'new'
			access = 'direct'
			form = 'unformatted'
			action = 'write'
			position = 'rewind'

		elseif(FileAccessType == PGSd_IO_Gen_USeqFrm) then
			status = 'old'
			access = 'sequential'
			form = 'formatted'
			action = 'readwrite'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_USeqUnf) then
			status = 'old'
			access = 'sequential'
			form = 'unformatted'
			action = 'readwrite'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_UDirFrm) then
			status = 'old'
			access = 'direct'
			form = 'formatted'
			action = 'readwrite'
			position = 'rewind'
		elseif(FileAccessType == PGSd_IO_Gen_UDirUnf) then
			status = 'old'
			access = 'direct'
			form = 'unformatted'
			action = 'readwrite'
			position = 'rewind'

		elseif(FileAccessType == PGSd_IO_Gen_ASeqFrm) then
			status = 'old'
			access = 'sequential'
			form = 'formatted'
			action = 'readwrite'
			position = 'append'
		elseif(FileAccessType == PGSd_IO_Gen_ASeqUnf) then
			status = 'old'
			access = 'sequential'
			form = 'unformatted'
			action = 'readwrite'
			position = 'append'

		else
			ErrType = UNKNOWNFILEACCESSTYPE
			return
		endif
		
		tiedup = .TRUE.

		do unit = bottom_unit_num, top_unit_num
      	inquire ( unit=unit, opened=tiedup )
			if (.not. tiedup) then
				exit
			endif
		enddo

		if(tiedup) THEN
			ErrType = NOFREEUNITS
			return
		endif
			
		print*, 'Fortran opening unit ', unit
		print*, 'access ', access
		print*, 'action ', action
		print*, 'form ', form
		print*, 'position ', position
		print*, 'status ', status

		if(access /= 'direct') then
			open(unit=unit, access=access, action=action, form=form, &
			& position=position, status=status, file=myName, iostat=ErrType)
		else
			open(unit=unit, access=access, action=action, form=form, &
			& status=status, file=myName, iostat=ErrType)
		endif
		
		print*, 'iostat ', ErrType
		call io_error('io error in MLSFiles: mls_io_gen_openF' // &
		& ' Fortran open', ErrType, myName)

		theFileHandle = unit
			
	case default
		ErrType = UNKNOWNTOOLBOXMODE

	end select
	
	if(ErrType /= 0) then
		theFileHandle = FH_ON_ERROR
	endif

  END FUNCTION mls_io_gen_openF

!====================
END MODULE MLSFiles
!====================

!
! $Log$
! Revision 2.4  2001/03/22 01:09:31  pwagner
! Added file name to Fortran open statement
!
! Revision 2.3  2001/03/21 00:48:43  pwagner
! Corrected mls_io_gen_openF
!
! Revision 2.2  2001/03/20 00:41:28  pwagner
! Added mls_io_gen_openF
!
! Revision 2.1  2001/03/07 01:02:37  pwagner
! First commit
!
!
