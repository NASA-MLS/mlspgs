
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============
MODULE Synoptic
!===============

   Use OpenInit
   Use MLSCommon
   USE MLSCF
   USE MLSPCF3
   USE L3CF
   USE L3DMData
   USE L2GPData
   USE L2Interface
   USE OutputClose
   Use global_data
   Use DailyMapModule
   USE L3DMData

   USE NumRecipesModule

   Implicit none

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- DailyCoreProcessing
!                
! This subroutine will process one product for all the input days data
!
! Remarks:  This module contains subroutines related to MLS L3 Daily Map Processing 
!

CONTAINS

!-------------------------------------------------------------------------
   SUBROUTINE DailyCoreProcessing(cfDef, cfProd, pcf, l2Days, l2gp, avgPeriod, l3sp, &
!     l3dm, dmA, dmD, l3r, residA, residD, flags)
     l3dm, dmA, dmD, l3r, residA, residD, mis_l2Days, mis_Days, flags)
!-------------------------------------------------------------------------

! Brief description of program
! This is the main program to run the Core processing.

! Parameters

	! Variable definitions

        TYPE( L3CFDef_T ) :: cfDef
        TYPE( PCFData_T ) :: pcf
        TYPE( L3CFProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
        TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
        TYPE( L3SPData_T ), POINTER :: l3sp(:), l3spPrec(:)
        TYPE( L2GPData_T ), POINTER :: l3r(:), residA(:), residD(:)

   	TYPE( OutputFlags_T ) :: flags

	Real (r8), POINTER, DIMENSION(:, :, :) ::  alons(:, :, :), alats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  dlons(:, :, :), dlats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  atimes(:, :, :), dtimes(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  aprec(:, :, :), dprec(:, :, :)

        REAL(r8), POINTER :: avgPeriod(:)

        INTEGER, DIMENSION(l2gp(1)%nLevels) :: perMisPoints
	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 
	Real (r8), DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: delTad 
        REAL(r8), DIMENSION(:), POINTER :: l3Result(:)       ! returned reconstructed result 

        REAL(r8), DIMENSION(:), POINTER :: sortTemp(:)       ! hold sort array temporary
        INTEGER, POINTER, DIMENSION(:) :: pt(:)

        INTEGER, POINTER, DIMENSION(:) :: nc(:), nca(:), ncd(:)
        Real (r8), POINTER, DIMENSION(:) :: startTime(:), endTime(:)

	INTEGER ::  error, l2Days, nlev, nlev_temp, nf, nwv, numDays, numSwaths
        INTEGER ::  rDays, pEndIndex, pStartIndex

        integer mis_l2Days, mis_l2Days_temp
        CHARACTER (LEN=DATE_LEN) :: mis_Days(maxWindow), mis_Days_temp(maxWindow)
        CHARACTER (LEN=480) :: msr

        integer i, j, iP, kP, iD, iL, iT
	real tau0, l3ret, avg
	real(r8) lonD0_in, lonA0_in


!*** Initilize variables
 
        nwv = cfProd%rangWavenumber(2) - cfProd%rangWavenumber(1) + 1
        nf = 60 

	nlev = 0
        DO j = 1, l2gp(1)%nLevels
           IF( l2gp(1)%pressures(j) >= cfProd%l3presLvl(1) .AND. &
	       l2gp(1)%pressures(j) <= cfProd%l3presLvl(2) ) THEN
	     nlev = nlev + 1
	     IF (nlev == 1) pStartIndex = j 
	     pEndIndex = j 
	   ENDIF
        ENDDO

	IF (nlev == 0) THEN
		nlev_temp = 1
		nlev = 1
                pStartIndex = 1
                pEndIndex = 1
	ELSE
		nlev_temp = -1
	END IF

!*** Initilize POINTERS

        IF (cfProd%mode == 'all' .or. cfProd%mode == 'ado') THEN
           numSwaths = 3
        ELSE
           numSwaths = 1
        ENDIF

        ALLOCATE( l3sp(numSwaths), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3sp array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( l3spPrec(numSwaths), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3spPrec array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        IF (cfProd%mode == 'all' .or. cfProd%mode == 'ado') THEN
           l3sp(1)%name = cfProd%l3prodNameD
           l3sp(2)%name = TRIM(cfProd%l3prodNameD) // 'Ascending'
           l3sp(3)%name = TRIM(cfProd%l3prodNameD) // 'Descending'
        ELSE
           l3sp(1)%name = TRIM(cfProd%l3prodNameD)
        ENDIF

        l3sp%startTime = l2gp(1)%time(1)
        l3sp%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)

        DO j = 1, numSwaths

           CALL AllocateL3SP( nlev, cfProd%nLats, nwv, nf, l3sp(j) )
           CALL AllocateL3SP( nlev, cfProd%nLats, nwv, nf, l3spPrec(j) )
           l3sp(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           l3sp(j)%latitude = cfProd%latGridMap(:l3sp(j)%nLats)
           l3sp(j)%waveNumber = 0.0
           l3sp(j)%frequency = 0.0
           l3sp(j)%l3spRelValue = 0.0 
           l3sp(j)%l3spRelPrecision = 0.0 
           l3sp(j)%l3spImgValue = 0.0 
           l3sp(j)%l3spImgPrecision = 0.0 

           l3spPrec(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           l3spPrec(j)%latitude = cfProd%latGridMap(:l3sp(j)%nLats)
           l3spPrec(j)%waveNumber = 0.0
           l3spPrec(j)%frequency = 0.0
           l3spPrec(j)%l3spRelValue = 0.0
           l3spPrec(j)%l3spRelPrecision = 0.0
           l3spPrec(j)%l3spImgValue = 0.0
           l3spPrec(j)%l3spImgPrecision = 0.0

        ENDDO

        numDays = cfProd%nDays

        if (numDays .gt. 0) then

        ALLOCATE( l3dm(numDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3dm array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( dmA(numDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' dmA array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( dmD(numDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' dmD array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        endif

!!      Initialize Daily Map & Diagnostic

        l3dm%name = cfProd%l3prodNameD
        dmA%name = TRIM(cfProd%l3prodNameD) // 'Ascending'
        dmD%name = TRIM(cfProd%l3prodNameD) // 'Descending'

        DO j = 1, numDays

           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, l3dm(j) )
           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmA(j) )
           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, cfDef%N, dmD(j) )

           l3dm(j)%time = cfProd%timeD(j)
           dmA(j)%time = cfProd%timeD(j)
           dmD(j)%time = cfProd%timeD(j)

           l3dm(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           dmA(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           dmD(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)

           l3dm(j)%latitude = cfProd%latGridMap(:l3dm(j)%nLats)
           dmA(j)%latitude = cfProd%latGridMap(:dmA(j)%nLats)
           dmD(j)%latitude = cfProd%latGridMap(:dmD(j)%nLats)

           l3dm(j)%longitude = cfProd%longGrid(:l3dm(j)%nLons)
           dmA(j)%longitude = cfProd%longGrid(:dmA(j)%nLons)
           dmD(j)%longitude = cfProd%longGrid(:dmD(j)%nLons)

           l3dm(j)%l3dmValue = 0.0
           dmA(j)%l3dmValue = 0.0
           dmD(j)%l3dmValue = 0.0

           l3dm(j)%l3dmPrecision = 0.0
           dmA(j)%l3dmPrecision = 0.0
           dmD(j)%l3dmPrecision = 0.0

           l3dm(j)%maxDiff = 0.0 
           dmA(j)%maxDiff = 0.0 
           dmD(j)%maxDiff = 0.0 

           l3dm(j)%maxDiffTime = 0.0 
           dmA(j)%maxDiffTime = 0.0 
           dmD(j)%maxDiffTime = 0.0 

        ENDDO

!!      Initialize Daily Map Residues

        CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, pcf%l3StartDay, & 
!          pcf%l3EndDay, rDays, l3r)
          pcf%l3EndDay, rDays, mis_l2Days_temp, mis_Days_temp, l3r)
        CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, pcf%l3StartDay, & 
!          pcf%l3EndDay, rDays, residA)
          pcf%l3EndDay, rDays, mis_l2Days_temp, mis_Days_temp, residA)
        CALL ReadL2GPProd(cfProd%l3prodNameD, cfProd%fileTemplate, pcf%l3StartDay, & 
!          pcf%l3EndDay, rDays, residD)
          pcf%l3EndDay, rDays, mis_l2Days_temp, mis_Days_temp, residD)

        l3r%name    = TRIM(cfProd%l3prodNameD) // 'Residuals'
        residA%name = TRIM(cfProd%l3prodNameD) // 'AscendingResiduals'
        residD%name = TRIM(cfProd%l3prodNameD) // 'DescendingResiduals'

        DO j = 1, rDays
           l3r(j)%l2gpValue    = 0.0
           l3r(j)%latitude     = 0.0
           l3r(j)%longitude    = 0.0
           residA(j)%l2gpValue = 0.0
           residA(j)%latitude  = 0.0
           residA(j)%longitude = 0.0
           residD(j)%l2gpValue = 0.0
           residD(j)%latitude  = 0.0
           residD(j)%longitude = 0.0
        ENDDO
 
!!      Initialize Flags

        flags%writel3sp    = .FALSE.
        flags%writel3dmCom = .FALSE.
        flags%writel3dmAsc = .FALSE.
        flags%writel3dmDes = .FALSE.
        flags%writel3rCom  = .FALSE.
        flags%writel3rAsc  = .FALSE.
        flags%writel3rDes  = .FALSE.

        if ( rDays .gt. 0) then 
 
        ALLOCATE( nc(rDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' nc array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF


        ALLOCATE( nca(rDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' nca array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( ncd(rDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' ncd array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( startTime(rDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' startTime array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        ALLOCATE( endTime(rDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' endTime array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        endif

        DO I = 1, rDays 
	    startTime(I) = l3r(I)%time(1)
	    endTime(I) = l3r(I)%time(l3r(I)%nTimes)
	ENDDO 
 
!!      Check if pressure levels are found 

	IF (nlev_temp == 1) THEN 
	   RETURN
	END IF

!*** Calculate average orbital period (day)

	tau0 = 0.0
        DO I = 1, size(avgPeriod)
	  tau0 = tau0 + avgPeriod(i)
	ENDDO 
        IF (size(avgPeriod) .ne. 0) tau0 = tau0/86400.0/float(size(avgPeriod))
 
!*** Sort & Prepare the Data 

	Call SortData(cfProd, l2Days, l2gp, 			&
			pStartIndex, pEndIndex,			&
			tau0, 					&
			anlats, dnlats, 			&
		        alats, dlats, 				&
			alons, dlons, 				&
			atimes, dtimes, 			&
			afields, dfields, 			&
			aprec, dprec,                           &
			delTad, perMisPoints )

!*** Main Loop 


   !*** Calculate Field Values ********

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1

          DO I = 1, rDays 
	    nc(I) = 0
	    nca(I) = 0
	    ncd(I) = 0
	  ENDDO 

	  DO J = 1, cfProd%nLats
       	     IF( anlats(J, iP) > 0 ) THEN

		IF( atimes(J, 1, iP) < dtimes(J, 1, iP) ) THEN

   		  lonD0_in = FindRealLon(real(dlons(J, 1, iP)))
   		  lonA0_in = FindRealLon(real(alons(J, 1, iP)))

	          CALL Init(cfProd%mode, 				&
			  nt_a_i   = anlats(J, iP), 			&
			  nt_d_i   = dnlats(J, iP), 			&
			  tau0_i   = tau0, 				&
			  delTad_i = delTad(J, iP), 			&
			  c0_i     = 2.0*PI, 				&
			  lonD0_i  = lonD0_in, 				&
			  tD0_i    = dtimes(J, 1, iP), 			&
			  lonA0_i  = lonA0_in, 				&
			  tA0_i    = atimes(J, 1, iP), 			&
			  lat_i    = alats(J, 1, iP) )
          		CALL DataGenerate(afields(J, :, iP), dfields(J, :, iP) )
		ELSE

   		  lonD0_in = FindRealLon(real(dlons(J, 2, iP)))
   		  lonA0_in = FindRealLon(real(alons(J, 1, iP)))

	          CALL Init(cfProd%mode, 				&
			  nt_a_i   = anlats(J, iP), 			&
			  nt_d_i   = dnlats(J, iP)-1, 			&
			  tau0_i   = tau0, 				&
			  delTad_i = delTad(J, iP), 			&
			  c0_i     = 2.0*PI, 				&
			  lonD0_i  = lonD0_in, 				&
			  tD0_i    = dtimes(J, 2, iP), 			&
			  lonA0_i  = lonA0_in, 				&
			  tA0_i    = atimes(J, 1, iP), 			&
			  lat_i    = alats(J, 1, iP) )
          		CALL DataGenerate(afields(J, :, iP), dfields(J, 2:, iP) )
		END IF

		IF (cfProd%mode == 'com') THEN

                  if (l3dm(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, &
                        real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, 	&
                        l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
          !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, l3dm(iD)%nLons
                        avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                    IF (l3dm(iD)%nLons .ne. 0) avg = avg/real(l3dm(iD)%nLons)
                      l3dm(iD)%latRss(iP, J) = 0.0
                      DO I = 1, l3dm(iD)%nLons
                        l3dm(iD)%latRss(iP, J) = l3dm(iD)%latRss(iP, J) +  &
                                      (l3dm(iD)%l3dmValue(iP, J, I)-avg)*  &
                                      (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO

                      IF (l3dm(iD)%nLons .ne. 0) THEN
                       l3dm(iD)%latRss(iP, J) = l3dm(iD)%latRss(iP, J)/ &
                         real(l3dm(iD)%nLons)
                      ENDIF

		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
		IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
                   atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
		   CALL Diagnostics(cfProd%mode, atimes(J, iL, iP), alons(J, iL, iP), l3ret) 
		   l3r(iD)%time(nc(iD))  = atimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
	           l3r(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
	           l3r(iD)%longitude(nc(iD)) = FindRealLon(real(alons(J, iL, iP)))*180.0/PI
		   l3r(iD)%l2gpValue(1, iP, nc(iD)) = afields(J, iL, iP)-l3ret 
	           l3r(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
                      DO iL = 1, dnlats(J, iP)
			IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
	       CALL Diagnostics(cfProd%mode, dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
		   l3r(iD)%time(nc(iD))     = dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
		   l3r(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
		   l3r(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
		   l3r(iD)%l2gpValue(1, iP, nc(iD)) = dfields(J, iL, iP)-l3ret 
		   l3r(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO
                   
                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmCom = .TRUE.
		   flags%writel3rCom  = .TRUE.
		   flags%writel3sp = .TRUE.

                   endif

		ELSE IF (cfProd%mode == 'asc') THEN

                  if (dmA(1)%nLons .gt. 0) then 

                   if ( associated(l3Result) ) then 
                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
       	      CALL Reconstruct(cfProd%mode, real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) +  &
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN 
                      dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons)
                      ENDIF
		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
			IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
		   CALL Diagnostics(cfProd%mode, atimes(J, iL, iP), alons(J, iL, iP), l3ret) 
	         residA(iD)%time(nc(iD))     = atimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
	         residA(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
	         residA(iD)%longitude(nc(iD)) = FindRealLon(real(alons(J, iL, iP)))*180.0/PI
	         residA(iD)%l2gpValue(1, iP, nc(iD)) = afields(J, iL, iP)-l3ret 
		 residA(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmAsc = .TRUE.
		   flags%writel3rAsc  = .TRUE.
		   flags%writel3sp = .TRUE.

                   endif

		ELSE IF (cfProd%mode == 'des') THEN

                  if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
             CALL Reconstruct(cfProd%mode, real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
			      dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) +  &
                                     (dmD(iD)%l3dmValue(iP, J, I)-avg)*  &
                                     (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN
                      dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J)/real(dmD(iD)%nLons)
                      ENDIF
		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
			IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
		   CALL Diagnostics(cfProd%mode, dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
		  residD(iD)%time(nc(iD))     = dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
		  residD(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
	          residD(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
	          residD(iD)%l2gpValue(1, iP, nc(iD)) = dfields(J, iL, iP)-l3ret 
		  residD(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmDes = .TRUE.
		   flags%writel3rDes  = .TRUE.
		   flags%writel3sp = .TRUE.

                   endif

                ELSE IF (cfProd%mode == 'ado') THEN

                  if (dmA(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform('asc')
                   CALL FFSMA(l3sp(2), iP, J)
                   DO iD = 1, cfProd%nDays
                   CALL Reconstruct('asc', real(dmA(iD)%time-l2gp(1)%time(1))/86400.0,  &
                                    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%l3dmValue(iP, J, I) = l3Result(I)
                      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) +  &
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)* & 
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN
                      dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons)
                      ENDIF
                   ENDDO 

                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
                        IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
                           atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
                           nca(iD) = nca(iD) + 1
                          CALL Diagnostics('asc', atimes(J, iL, iP), alons(J, iL, iP), l3ret)
                           residA(iD)%time(nc(iD))     = atimes(J, iL, iP)*86400.0+ & 
                             l2gp(1)%time(1)
                           residA(iD)%latitude(nca(iD)) = cfProd%latGridMap(J)  
                           residA(iD)%longitude(nc(iD)) = &
                             FindRealLon(real(alons(J, iL, iP)))*180.0/PI
                           residA(iD)%l2gpValue(1, iP, nca(iD)) = afields(J, iL, iP)-l3ret
                           residA(iD)%l2gpPrecision(1, iP, nca(iD)) = 0.0
                        END IF
                      ENDDO
                   ENDDO

                   if ( associated(l3Result) ) then 
                   DeAllocate(l3Result, STAT=error) 
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   flags%writel3dmAsc = .TRUE. 
                   flags%writel3rAsc  = .TRUE.

                   endif

                   if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform('des')
                   CALL FFSMD(l3sp(3), iP, J)
                   DO iD = 1, cfProd%nDays
                      CALL Reconstruct('des', real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
                                    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%l3dmValue(iP, J, I) = l3Result(I)
                      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) + &
                                               (dmD(iD)%l3dmValue(iP, J, I)-avg)* & 
                                               (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN 
                      dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J)/ &
                        real(dmD(iD)%nLons)
                      ENDIF
                   ENDDO 

                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
                        IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
                           dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
                           ncd(iD) = ncd(iD) + 1
                          CALL Diagnostics('des', dtimes(J, iL, iP), dlons(J, iL, iP), l3ret)
                  residD(iD)%time(ncd(iD))     = dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
                  residD(iD)%latitude(ncd(iD)) = cfProd%latGridMap(J)  
                  residD(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
                  residD(iD)%l2gpValue(1, iP, ncd(iD)) = dfields(J, iL, iP)-l3ret
                  residD(iD)%l2gpPrecision(1, iP, ncd(iD)) = 0.0
                        END IF
                      ENDDO
                   ENDDO

                   if ( associated(l3Result) ) then 
                   DeAllocate(l3Result, STAT=error) 
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   flags%writel3dmDes = .TRUE. 
                   flags%writel3rDes  = .TRUE.
                   flags%writel3sp = .TRUE.

                   endif

		ELSE IF (cfProd%mode == 'all') THEN

                  if (l3dm(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform('com')
	   	   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
            CALL Reconstruct('com', real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, &
				    l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, l3dm(iD)%nLons
                        avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (l3dm(iD)%nLons .ne. 0) avg = avg/real(l3dm(iD)%nLons)
                      l3dm(iD)%latRss(iP, J) = 0.0
                      DO I = 1, l3dm(iD)%nLons
                        l3dm(iD)%latRss(iP, J) = l3dm(iD)%latRss(iP, J) +  &
                                      (l3dm(iD)%l3dmValue(iP, J, I)-avg)*  &
                                      (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (l3dm(iD)%nLons .ne. 0) THEN
                      l3dm(iD)%latRss(iP, J) = l3dm(iD)%latRss(iP, J)/ &
                        real(l3dm(iD)%nLons)
                      ENDIF
		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
			IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
		         CALL Diagnostics('com', atimes(J, iL, iP), alons(J, iL, iP), l3ret) 
		        l3r(iD)%time(nc(iD))  = atimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
			l3r(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
			l3r(iD)%longitude(nc(iD)) = & 
                          FindRealLon(real(alons(J, iL, iP)))*180.0/PI
			   l3r(iD)%l2gpValue(1, iP, nc(iD)) = afields(J, iL, iP)-l3ret 
			   l3r(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
                      DO iL = 1, dnlats(J, iP)
			IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nc(iD) = nc(iD) + 1
		         CALL Diagnostics('com', dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
	             l3r(iD)%time(nc(iD))     = dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
	             l3r(iD)%latitude(nc(iD)) = cfProd%latGridMap(J) 
		     l3r(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
		     l3r(iD)%l2gpValue(1, iP, nc(iD)) = dfields(J, iL, iP)-l3ret 
		     l3r(iD)%l2gpPrecision(1, iP, nc(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result,STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmCom = .TRUE.
		   flags%writel3rCom  = .TRUE.

                   endif

                   if (dmA(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3sp(2), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('asc', real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, &
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) avg = avg/real(dmA(iD)%nLons)
                      dmA(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J) + &
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                                               (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmA(iD)%nLons .ne. 0) THEN
                      dmA(iD)%latRss(iP, J) = dmA(iD)%latRss(iP, J)/real(dmA(iD)%nLons)
                      ENDIF
		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, anlats(J, iP)
			IF(atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   atimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   nca(iD) = nca(iD) + 1
	           CALL Diagnostics('asc', atimes(J, iL, iP), alons(J, iL, iP), l3ret) 
		  residA(iD)%time(nc(iD))     = atimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
		  residA(iD)%latitude(nca(iD)) = cfProd%latGridMap(J) 
		  residA(iD)%longitude(nc(iD)) = FindRealLon(real(alons(J, iL, iP)))*180.0/PI
                  residA(iD)%l2gpValue(1, iP, nca(iD)) = afields(J, iL, iP)-l3ret 
                  residA(iD)%l2gpPrecision(1, iP, nca(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result,STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmAsc = .TRUE.
		   flags%writel3rAsc  = .TRUE.

                   endif

                   if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3sp(3), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('des', real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, &
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		      ENDDO
                      !***Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)
                      avg = 0.0
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) avg = avg/real(dmD(iD)%nLons)
                      dmD(iD)%latRss(iP, J) = 0.0
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J) +  &
                                               (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                                               (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      IF (dmD(iD)%nLons .ne. 0) THEN 
                      dmD(iD)%latRss(iP, J) = dmD(iD)%latRss(iP, J)/real(dmD(iD)%nLons)
                      ENDIF
		   ENDDO

                   DO iD = 1, rDays
                      DO iL = 1, dnlats(J, iP)
			IF(dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) >= startTime(iD) .AND. &
			   dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1) <= endTime(iD)) THEN
			   ncd(iD) = ncd(iD) + 1
		         CALL Diagnostics('des', dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
		  residD(iD)%time(ncd(iD))     = dtimes(J, iL, iP)*86400.0+l2gp(1)%time(1)
                  residD(iD)%latitude(ncd(iD)) = cfProd%latGridMap(J) 
                  residD(iD)%longitude(nc(iD)) = FindRealLon(real(dlons(J, iL, iP)))*180.0/PI
                  residD(iD)%l2gpValue(1, iP, ncd(iD)) = dfields(J, iL, iP)-l3ret 
                  residD(iD)%l2gpPrecision(1, iP, ncd(iD)) = 0.0 
			END IF
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result,STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

		   flags%writel3dmDes = .TRUE.
		   flags%writel3rDes  = .TRUE.
		   flags%writel3sp = .TRUE.

                   endif

		END IF

                CALL ClearMemory()

	     END IF
	  ENDDO

          !*** Sort into time ascending order according to l2gp format
          IF (cfProd%mode == 'com') THEN
                DO iD = 1, rDays

                  if (nc(iD) .gt. 0) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      DO i = 1, nc(iD)
                        pt(i) = 0
                        sortTemp(i) = 0.0
                      ENDDO
                      CALL DSORTP(l3r(iD)%time, 1, nc(iD), pt)
                      !** do time
                      CALL DSORT(l3r(iD)%time, 1, nc(iD))
                      !** do latitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO
                      
                      if ( associated(sortTemp) ) then 
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      if ( associated(pt)) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if ( nc(iD) .gt. 0 ) then 
                      
                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(l3r(iD)%l2gpValue(1, iP, 1:nc(iD))), 1, nc(iD), pt)
                      DO i = 1, cfDef%N
                        l3dm(iD)%maxDiff(i, iP) = l3r(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, l3dm(iD)%nLons
                        avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      avg = avg/real(l3dm(iD)%nLons*cfProd%nLats)
                      ENDIF
                      l3dm(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, l3dm(iD)%nLons
                        l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP) +         &
                                           (l3dm(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP)/ &
                        real(l3dm(iD)%nLons*cfProd%nLats)
                      ENDIF
                      l3dm(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
          ELSE IF (cfProd%mode == 'asc') THEN
                DO iD = 1, rDays

                  if ( nc(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residA(iD)%time, 1, nc(iD), pt)
                      !** do time
                      CALL DSORT(residA(iD)%time, 1, nc(iD))
                      !** do latitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = residA(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residA(iD)%latitude(i) = sortTemp(pt(i))
                      !** do longitude
                      ENDDO
                      DO i = 1, nc(iD)
                        sortTemp(i) = residA(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residA(iD)%longitude(i) = sortTemp(pt(i))
                      !** do value
                      ENDDO
                      DO i = 1, nc(iD)
                        sortTemp(i) = residA(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residA(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO

                      if ( associated(sortTemp)) then
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      if ( associated(pt) ) then
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if ( nc(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residA(iD)%l2gpValue(1, iP, :)), 1, nc(iD), pt)
                      DO i = 1,cfDef%N
                        dmA(iD)%maxDiff(i, iP) = residA(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmA(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP)/ &
                        real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
          ELSE IF (cfProd%mode == 'des') THEN
                DO iD = 1, rDays

                  if ( nc(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residD(iD)%time, 1, nc(iD), pt)
                      !** do time
                      CALL DSORT(residD(iD)%time, 1, nc(iD))
                      !** do latitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = residD(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residD(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = residD(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residD(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nc(iD)
                        sortTemp(i) = residD(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nc(iD)
                        residD(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO

                      if ( associated(sortTemp) ) then 
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if (nc(iD) .gt. 0) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residD(iD)%l2gpValue(1, iP, :)), 1, nc(iD), pt)
                      DO i = 1,cfDef%N
                        dmD(iD)%maxDiff(i, iP) = residD(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP)/ & 
                        real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
          ELSE IF (cfProd%mode == 'ado') THEN
                DO iD = 1, rDays

                  if ( nca(iD) .gt. 0 ) then 
                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residA(iD)%time, 1, nca(iD), pt)
                      !** do time
                      CALL DSORT(residA(iD)%time, 1, nca(iD))
                      !** do latitude
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO
                      if ( associated(sortTemp) ) then 
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif
                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference
                      
                      if ( nca(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residA(iD)%l2gpValue(1, iP, :)), 1, nca(iD), pt)
                      DO i = 1,cfDef%N
                        dmA(iD)%maxDiff(i, iP) = residA(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO

                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                        avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF

                      dmA(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO

                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0 ) THEN 
                      dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP)/ & 
                        real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF

                      dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO

                DO iD = 1, rDays

                  if (ncd(iD) .gt. 0) then 
                  
                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residD(iD)%time, 1, ncd(iD), pt)
                      !** do time
                      CALL DSORT(residD(iD)%time, 1, ncd(iD))
                      !** do latitude
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%latitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%longitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO

                      if ( associated(sortTemp) ) then 
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif
                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if (ncd(iD) .gt. 0) then 

                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residD(iD)%l2gpValue(1, iP, :)), 1, ncd(iD), pt)
                      DO i = 1,cfDef%N
                        dmD(iD)%maxDiff(i, iP) = residD(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN
                      avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP)/ & 
                        real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
          ELSE IF (cfProd%mode == 'all') THEN
                DO iD = 1, rDays

                  if (nca(iD) .gt. 0) then 

                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residA(iD)%time, 1, nca(iD), pt)
                      !** do time
                      CALL DSORT(residA(iD)%time, 1, nca(iD))
                      !** do latitude
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nca(iD)
                        sortTemp(i) = residA(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nca(iD)
                        residA(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO

                      if ( associated(sortTemp)) then 

                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      endif

                      if ( associated(pt) ) then 

                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if (nca(iD) .gt. 0) then 

                      ALLOCATE( pt(nca(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residA(iD)%l2gpValue(1, iP, :)), 1, nca(iD), pt)
                      DO i = 1,cfDef%N
                        dmA(iD)%maxDiff(i, iP) = residA(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        avg = avg + dmA(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      avg = avg/real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmA(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmA(iD)%nLons
                        dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP) +   &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmA(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (dmA(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmA(iD)%gRss(iP) = dmA(iD)%gRss(iP)/ & 
                        real(dmA(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmA(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO

                DO iD = 1, rDays

                  if (ncd(iD) .gt. 0) then 

                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(residD(iD)%time, 1, ncd(iD), pt)
                      !** do time
                      CALL DSORT(residD(iD)%time, 1, ncd(iD))
                      !** do latitude
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%latitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%longitude(i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, ncd(iD)
                        sortTemp(i) = residD(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, ncd(iD)
                        residD(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO
                      if ( associated(sortTemp) ) then
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif
                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if ( ncd(iD) .gt. 0 ) then 

                      ALLOCATE( pt(ncd(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(residD(iD)%l2gpValue(1, iP, :)), 1, ncd(iD), pt)
                      DO i = 1,cfDef%N
                        dmD(iD)%maxDiff(i, iP) = residD(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        avg = avg + dmD(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      avg = avg/real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, dmD(iD)%nLons
                        dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP) +   &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (dmD(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (dmD(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      dmD(iD)%gRss(iP) = dmD(iD)%gRss(iP)/real(dmD(iD)%nLons*cfProd%nLats)
                      ENDIF
                      dmD(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO

                DO iD = 1, rDays

                  if ( nc(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      ALLOCATE( sortTemp(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(l3r(iD)%time, 1, nc(iD), pt)
                      !** do time
                      CALL DSORT(l3r(iD)%time, 1, nc(iD))
                      !** do latitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%latitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%latitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do longitude
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%longitude(i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%longitude(i) = sortTemp(pt(i))
                      ENDDO
                      !** do value
                      DO i = 1, nc(iD)
                        sortTemp(i) = l3r(iD)%l2gpValue(1, iP, i)
                      ENDDO
                      DO i = 1, nc(iD)
                        l3r(iD)%l2gpValue(1, iP, i) = sortTemp(pt(i))
                      ENDDO

                      if ( associated(sortTemp) ) then 
                      DeAllocate(sortTemp, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  sortTemp array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !*** Calculate Maxmum Difference

                      if ( nc(iD) .gt. 0 ) then 

                      ALLOCATE( pt(nc(iD)), STAT=error )
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_Allocate // ' pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF

                      CALL DSORTP(-abs(l3r(iD)%l2gpValue(1, iP, :)), 1, nc(iD), pt)
                      DO i = 1,cfDef%N
                        l3dm(iD)%maxDiff(i, iP) = l3r(iD)%l2gpValue(1, iP, pt(i))
                      ENDDO

                      if ( associated(pt) ) then 
                      DeAllocate(pt, STAT=error)
                      IF ( error /= 0 ) THEN
                        msr = MLSMSG_DeAllocate // '  pt array.'
                        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                      ENDIF
                      endif

                      endif

                      !***Root-Sum-Square for each latitude, dimensioned (nLevels)
                      avg = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, l3dm(iD)%nLons
                        avg = avg + l3dm(iD)%l3dmValue(iP, J, I)
                      ENDDO
                      ENDDO
                      IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                        avg = avg/real(l3dm(iD)%nLons*cfProd%nLats)
                      ENDIF
                      l3dm(iD)%gRss(iP) = 0.0
                      DO J = 1, cfProd%nLats
                      DO I = 1, l3dm(iD)%nLons
                        l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP) +         &
                                           (l3dm(iD)%l3dmValue(iP, J, I)-avg)* &
                                           (l3dm(iD)%l3dmValue(iP, J, I)-avg)
                      ENDDO
                      ENDDO
                      IF (l3dm(iD)%nLons*cfProd%nLats .ne. 0) THEN 
                      l3dm(iD)%gRss(iP) = l3dm(iD)%gRss(iP)/ &
                        real(l3dm(iD)%nLons*cfProd%nLats)
                      ENDIF
                      l3dm(iD)%perMisPoints(iP) = perMisPoints(iP)
                ENDDO
          END IF

        ENDDO

   !*** Calculate Field Precisions ********

	iP = 0
	DO kP = pStartIndex, pStartIndex 
	  iP = iP + 1

	  DO J = 1, cfProd%nLats
       	     IF( anlats(J, iP) > 0 ) THEN

		IF( atimes(J, 1, iP) < dtimes(J, 1, iP) ) THEN

   		  lonD0_in = FindRealLon(real(dlons(J, 1, iP)))
   		  lonA0_in = FindRealLon(real(alons(J, 1, iP)))

	          CALL Init(cfProd%mode, 					&
			  nt_a_i   = anlats(J, iP), 				&
			  nt_d_i   = dnlats(J, iP), 				&
			  tau0_i   = tau0, 					&
			  delTad_i = delTad(J, iP), 				&
			  c0_i     = 2.0*PI, 					&
			  lonD0_i  = lonD0_in, 					&
			  tD0_i    = dtimes(J, 1, iP), 				&
			  lonA0_i  = lonA0_in, 					&
			  tA0_i    = atimes(J, 1, iP), 				&
			  lat_i    = alats(J, 1, iP) )
          	  CALL DataGenerate(aprec(J, :, iP), dprec(J, :, iP) )

		ELSE

   		  lonD0_in = FindRealLon(real(dlons(J, 2, iP)))
   		  lonA0_in = FindRealLon(real(alons(J, 1, iP)))

	          CALL Init(cfProd%mode, 					&
			  nt_a_i   = anlats(J, iP), 				&
			  nt_d_i   = dnlats(J, iP)-1, 				&
			  tau0_i   = tau0, 					&
			  delTad_i = delTad(J, iP), 				&
			  c0_i     = 2.0*PI, 					&
			  lonD0_i  = lonD0_in, 					&
			  tD0_i    = dtimes(J, 2, iP), 				&
			  lonA0_i  = lonA0_in, 					&
			  tA0_i    = atimes(J, 1, iP), 				&
			  lat_i    = alats(J, 1, iP) )
         	  CALL DataGenerate(aprec(J, :, iP), dprec(J, 2:, iP) )

		END IF

		IF (cfProd%mode == 'com') THEN

                  if (l3dm(1)%nLons .gt. 0) then 
                  
                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, &
                        real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

		ELSE IF (cfProd%mode == 'asc') THEN

                  if (dmA(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, &
                        real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

		ELSE IF (cfProd%mode == 'des') THEN

                  if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, &
                        real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

		ELSE IF (cfProd%mode == 'ado') THEN

                  if (dmA(1)%nLons .gt. 0) then 
                  
                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3spPrec(2), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('asc', &
                        real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

                   if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3spPrec(3), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('des', &
                        real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

		ELSE IF (cfProd%mode == 'all') THEN

                  if (l3dm(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform('com')
	   	   CALL FFSM(l3spPrec(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('com', &
                        real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

                   if (dmA(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF

                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3spPrec(2), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('asc', &
                        real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

                   if (dmD(1)%nLons .gt. 0) then 

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_Allocate // ' l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3spPrec(3), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('des', &
                        real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmPrecision(iP, J, I) = l3Result(I) 
		      ENDDO
		   ENDDO

                   if ( associated(l3Result) ) then 
		   DeAllocate(l3Result, STAT=error)
                   IF ( error /= 0 ) THEN
                     msr = MLSMSG_DeAllocate // '  l3Result array.'
                     CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
                   ENDIF
                   endif

                   endif

		END IF

                CALL ClearMemory()

	     END IF
	  ENDDO

	ENDDO

        if ( associated(l3spPrec) ) then 
	DeAllocate(l3spPrec, STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_DeAllocate // '  l3spPrec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(nc) ) then
	DeAllocate(nc, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  nc array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(nca) ) then 
	DeAllocate(nca, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  nca array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(ncd) ) then
	DeAllocate(ncd, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  ncd array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(startTime)) then 
	DeAllocate(startTime, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  startTime array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(endTime)) then
	DeAllocate(endTime, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  endTime array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(alats) ) then
	DeAllocate(alats, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  alats array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(dlats) ) then 
	DeAllocate(dlats, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  dlats array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(alons)) then
	DeAllocate(alons, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  alons array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(dlons)) then
	DeAllocate(dlons, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  dlons array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(atimes)) then
	DeAllocate(atimes, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  atimes array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(dtimes) ) then
	DeAllocate(dtimes, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  dtimes array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(afields) ) then
	DeAllocate(afields, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  afields array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(dfields) ) then
	DeAllocate(dfields, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  dfields array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(aprec)) then 
	DeAllocate(aprec, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  aprec array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(dprec)) then
	DeAllocate(dprec, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  dprec array.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

!-----------------------------------
   END SUBROUTINE DailyCoreProcessing
!-----------------------------------



!-------------------------------------------------------------------------
   SUBROUTINE SortData(cfProd, l2Days, l2gp, pStartIndex, pEndIndex, tau0, anlats, dnlats, &
		alats_interp, dlats_interp, alons_interp, dlons_interp, 	&
		atimes_interp, dtimes_interp, afields_interp, dfields_interp,	&
		aprec_interp, dprec_interp, delTad, perMisPoints)
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex_prev

        TYPE( L3CFProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	Real (r8), POINTER, DIMENSION(:, :) ::  &
          l2Times(:, :), l2Lons_new(:, :), &
          l2Lons(:, :),  l2Lons_old(:, :), l2Lats(:, :), l2Values(:, :), l2Prec(:, :)

	Real (r8), POINTER, DIMENSION(:, :, :) ::  & 
          alons_interp(:, :, :), alats_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  & 
          dlons_interp(:, :, :), dlats_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  & 
          atimes_interp(:, :, :), dtimes_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  & 
          afields_interp(:, :, :), dfields_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  &
          aprec_interp(:, :, :), dprec_interp(:, :, :)

	Real (r8) &
          dlons, alons, lons_found, lats_found, times_found, fields_found, prec_found, 	&
	  slope, slope_time, slope_field, slope_prec, sTime

	Real (r8) lons_found_old

	Real tau0
 
	Integer nPd, pStartIndex, pEndIndex, iMin

	INTEGER, DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: anlats, dnlats
	Real (r8), DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: delTad 

        INTEGER, DIMENSION(pEndIndex-pStartIndex+1) :: perMisPoints, numData
        CHARACTER (LEN=480) :: msr

!*** Calculate the number of points in each pressure level 

	nPd = 15

        nterms = 0
	DO I = 1, l2Days
            nterms = nterms + l2gp(I)%nTimes
	ENDDO

!*** Allocate space for all the arrays

        if ((nterms.gt.0).and.( (pEndIndex-pStartIndex).ge.0) ) then 

        ALLOCATE(l2Times(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Times array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF


        ALLOCATE(l2Lons(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lons array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF


        ALLOCATE(l2Lats(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lats array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF


	ALLOCATE(l2Values(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Values array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF


        ALLOCATE(l2Lons_old(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lons_old array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

	ALLOCATE(l2Prec(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Prec array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

        ALLOCATE(l2Lons_new(nterms, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' l2Lons_new array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          endif

        if ( (cfProd%nLats.gt.0).and.( (pEndIndex-pStartIndex).ge.0).and.(nPd*l2Days.gt.0)  ) then 

        ALLOCATE(alons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' alons_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

        ALLOCATE(alats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' alats_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

       ALLOCATE(atimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' atimes_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

      ALLOCATE(afields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' afields_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

	ALLOCATE(aprec_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' aprec_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

        ALLOCATE(dlons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dlons_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

        ALLOCATE(dlats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dlats_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

       ALLOCATE(dtimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dtimes_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

      ALLOCATE(dfields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dfields_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

	ALLOCATE(dprec_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

          IF ( error /= 0 ) THEN
             msr = MLSMSG_Allocate // ' dprec_interp array.'
             CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF

          endif


!*** Re-arrange the data into longitude order for each pressure level 

        DO kP = pStartIndex, pEndIndex
          perMisPoints(kP+1-pStartIndex) = 0
          numData(kP+1-pStartIndex) = 0
        ENDDO

        nstart = 1
	DO iD = 1, l2Days
	  DO iT = 1, l2gp(iD)%nTimes
	    iP = 0
	    DO kP = pStartIndex, pEndIndex 
	        iP = iP + 1
                l2Times(nstart, iP) = l2gp(iD)%time(iT)
                l2Lons(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lons_old(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lats(nstart, iP) = l2gp(iD)%latitude(iT)
                l2Values(nstart, iP) = l2gp(iD)%l2gpValue(1, kP, iT) 
		l2Prec(nstart, iP) = l2gp(iD)%l2gpPrecision(1, kP, iT)
                IF(l2Values(nstart, iP) < 0.0) perMisPoints(iP) = perMisPoints(iP) + 1
                numData(iP) = numData(iP) + 1
	    ENDDO
            nstart = nstart + 1
	  ENDDO
	ENDDO

        DO kP = pStartIndex, pEndIndex
          if (numData(kP+1-pStartIndex).ne.0) then 
            perMisPoints(kP+1-pStartIndex) = & 
              100*perMisPoints(kP+1-pStartIndex)/numData(kP+1-pStartIndex)
            endif
        ENDDO

!*** Re-arrange the longitude into sequence 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
          nr = 0
	  DO iT = 2, nterms 
                IF( l2Lons(iT, iP) >= 0.0 .AND. l2Lons(iT-1, iP) < 0.0 ) THEN 
     			l2Lons_new(iT, iP) = l2Lons(iT, iP) - nr*2.0*PI - 2.0*PI
     			nr = nr + 1
		ELSE
     			l2Lons_new(iT, iP) = l2Lons(iT, iP) - nr*2.0*PI
		END IF
	  ENDDO
          l2Lons_new(1, iP) = l2Lons(1, iP)
	ENDDO

!*** Convert time reference to the starting point

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
          sTime = l2Times(1, iP)
	  DO iT = 1, nterms 
                l2Times(iT, iP) = (l2Times(iT, iP) - sTime)/86400.0 
	  ENDDO
	ENDDO

!*** Find longitudes & latitudes for L3 latitudes in L2 data points 

        dlons = 0.3*PI/180.0

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1

 	  alons = l2Lons_new(1, iP) 

	  lindex = 1
	  lindex_prev = 0
	  aindex_prev = 0
	  dindex_prev = 0

          nloop = (l2Lons_new(1, iP)-l2Lons_new(nterms, iP))/dlons

	  DO J = 1, cfProd%nLats
       	     anlats(J, iP) = 0 
       	     dnlats(J, iP) = 0 
	  ENDDO

	  DO I = 1, nloop
	    alons = alons -  dlons
	    DO iT = lindex, nterms-1 
              IF( l2Lons_new(iT, iP) > alons .AND. l2Lons_new(iT+1, iP) <= alons ) THEN
	 	 lindex = iT
	 	 EXIT
	      END IF
	    ENDDO

            IF ( lindex /= lindex_prev ) THEN

	     lindex_prev = lindex

	     DO J = 1, cfProd%nLats
              IF( ( l2Lats(lindex, iP) <= cfProd%latGridMap(J) .AND. 	&
		    l2Lats(lindex+1, iP) > cfProd%latGridMap(J) )  .OR.	&
                  ( l2Lats(lindex, iP) >= cfProd%latGridMap(J) .AND. 	&
		    l2Lats(lindex+1, iP) < cfProd%latGridMap(J) )  ) THEN
                
                if (abs(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP)).ge.1.0e-12) then 
		slope = (l2Lats(lindex+1, iP)-l2Lats(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
                endif

                if (abs(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP)).ge.1.0e-12) then 
		slope_time = (l2Times(lindex+1, iP)-l2Times(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
                endif

                if (abs(l2Lons(lindex+1, iP)-l2Lons(lindex, iP)).ge.1.0e-12) then 
		slope_field = (l2Values(lindex+1, iP)-l2Values(lindex, iP))/	&
			(l2Lons(lindex+1, iP)-l2Lons(lindex, iP))
                endif 

                if (abs(l2Lons(lindex+1, iP)-l2Lons(lindex, iP)).ge.1.0e-12) then 
                slope_prec = (l2Prec(lindex+1, iP)-l2Prec(lindex, iP))/ &
                        (l2Lons(lindex+1, iP)-l2Lons(lindex, iP))
                endif

		IF(ABS(slope) < 1.e-6) THEN
		   lons_found = l2Lons_new(lindex, iP)
		   lons_found_old = l2Lons_old(lindex, iP)
		ELSE
		   lons_found = l2Lons_new(lindex, iP) + & 
                     (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		   lons_found_old = l2Lons_old(lindex, iP) + &
                     (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		END IF
		lats_found = cfProd%latGridMap(J)
		times_found = l2Times(lindex, iP) + &
                  (lons_found-l2Lons_new(lindex, iP))*slope_time
		fields_found = l2Values(lindex, iP) + & 
                  (lons_found-l2Lons_new(lindex, iP))*slope_field
		prec_found = l2Prec(lindex, iP) + & 
                  (lons_found-l2Lons_new(lindex, iP))*slope_prec
                aindex = J
		EXIT
	      END IF
	     ENDDO
              
	     IF( slope < 0 .and. aindex /= aindex_prev) THEN
	         aindex_prev = aindex
       	           anlats(aindex, iP) = anlats(aindex, iP) + 1 
	           alons_interp(aindex, anlats(aindex, iP), iP) = lons_found
	           alats_interp(aindex, anlats(aindex, iP), iP) = lats_found
	           atimes_interp(aindex, anlats(aindex, iP), iP) = times_found
	           afields_interp(aindex, anlats(aindex, iP), iP) = fields_found
		   aprec_interp(aindex, anlats(aindex, iP), iP) = prec_found
	     END IF

	     IF( slope >= 0 .and. aindex /= dindex_prev) THEN
	         dindex_prev = aindex
		 IF(anlats(aindex, iP) > 0) THEN
       	           dnlats(aindex, iP) = dnlats(aindex, iP) + 1 
	           dlons_interp(aindex, dnlats(aindex, iP), iP) = lons_found
	           dlats_interp(aindex, dnlats(aindex, iP), iP) = lats_found
	           dtimes_interp(aindex, dnlats(aindex, iP), iP) = times_found
	           dfields_interp(aindex, dnlats(aindex, iP), iP) = fields_found
		   dprec_interp(aindex, dnlats(aindex, iP), iP) = prec_found
	         END IF
	     END IF

	    END IF

          ENDDO

        ENDDO

!*** Find averaged asc/des time difference 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO J = 1, cfProd%nLats
            delTad(J,iP) = 0.0
            IF (anlats(J, iP) < dnlats(J, iP)) THEN
		iMin = anlats(J, iP)
	    ELSE
		iMin = dnlats(J, iP)
	    END IF 
	    DO I = 2, iMin 
             delTad(J,iP) = delTad(J,iP) + dtimes_interp(J, I, iP) - atimes_interp(J, I, iP)
	    ENDDO
            IF(iMin > 0) THEN
                if (iMin .ne. 1) delTad(J,iP) = delTad(J,iP)/float(iMin-1)
	    ELSE
                delTad(J,iP) = 0.0 
	    END IF

            IF(delTad(J,iP) < 0) THEN
                delTad(J,iP) = delTad(J,iP) + tau0 
	    END IF
	  ENDDO
	ENDDO

!*** Deallocate intermediate arrays 

        if ((nterms.gt.0).and.( (pEndIndex-pStartIndex).ge.0) ) then 

        if ( associated(l2Times) ) then
        DEALLOCATE(l2Times, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Times pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Lons) ) then 
        DEALLOCATE(l2Lons, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Lons pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Lons_new) ) then
        DEALLOCATE(l2Lons_new, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Lons_new pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Lons_old) ) then
        DEALLOCATE(l2Lons_old, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Lons_old pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Lats) ) then
        DEALLOCATE(l2Lats, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Lats pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Values) ) then
        DEALLOCATE(l2Values, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Values pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        if ( associated(l2Prec) ) then
        DEALLOCATE(l2Prec, STAT=error)
          IF ( error /= 0 ) THEN
              msr = MLSMSG_DeAllocate // '  l2Prec pointer.'
              CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          ENDIF
        endif

        endif

!-----------------------------------
   END SUBROUTINE SortData
!-----------------------------------


!-------------------------------------------------------------------------
   REAL(r8) FUNCTION FindRealLon(alon)
!-------------------------------------------------------------------------

	Real alon

	IF(mod(alon, 2.0*PI) < -PI) THEN
	   FindRealLon = 2.0*PI + mod(alon, 2.0*PI) 
	ELSE IF(mod(alon, 2.0*PI) > PI) THEN
	   FindRealLon = -2.0*PI + mod(alon, 2.0*PI) 
	ELSE 
	   FindRealLon = mod(alon, 2.0*PI) 
	END IF

!-----------------------------------
   END FUNCTION FindRealLon
!-----------------------------------



!===================
END MODULE Synoptic
!===================

! $Log$
! Revision 1.23  2002/04/01 21:57:13  jdone
! check division by zero
!
! Revision 1.22  2002/03/27 22:45:17  jdone
! deallocate statements added
!
! Revision 1.21  2002/03/27 21:35:14  jdone
! allocate statements checked and maxDiff is initialized
!
! Revision 1.20  2002/02/20 22:30:33  ybj
! *** empty log message ***
!
! Revision 1.19  2002/02/20 22:13:54  ybj
! *** empty log message ***
!
! Revision 1.18  2001/09/27 20:38:31  ybj
! Add Precision Calculation & ado mode
!
! Revision 1.17  2001/09/27 20:07:14  ybj
! Add Precision Calculation
!
! Revision 1.16  2001/09/06 18:45:37  nakamura
! Initialized pStart/EndIndex; removed some unused variables.
!
! Revision 1.15  2001/08/13 19:28:03  ybj
! *** empty log message ***
!
! Revision 1.14  2001/08/13 16:42:20  ybj
! *** empty log message ***
!
! Revision 1.13  2001/04/12 16:04:41  ybj
! reasonable values
!
! Revision 1.12  2001/04/11 18:29:22  ybj
! reasonable values
!
! Revision 1.11  2001/03/08 00:53:15  ybj
! *** empty log message ***
!
! Revision 1.10  2001/03/08 00:37:17  ybj
! *** empty log message ***
!
! Revision 1.9  2001/03/07 23:12:15  ybj
! *** empty log message ***
!
! Revision 1.8  2001/03/05 19:58:59  ybj
! with selected pressure levels
!
! Revision 1.7  2001/03/03 01:38:02  ybj
! with selected pressure levels
!
! Revision 1.6  2001/03/03 00:45:44  ybj
! with selected pressure levels
!
! Revision 1.5  2001/03/03 00:20:07  ybj
! with selected pressure levels
!
! Revision 1.3  2001/02/28 19:27:03  ybj
! Fix residue allocation & flags initialization
!
! Revision 1.2  2001/02/28 17:15:06  ybj
! Fix residue structure alocation
!
! Revision 1.1  2001/02/27 20:51:28  ybj
! Daily Map Core Processing
!
! Revision 1.1  2000/10/05 18:38:48  nakamura
! Program split from synoptic.f90 and modified to be more like the standard template.
!


