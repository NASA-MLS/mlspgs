
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============
MODULE MonthlyProcessModule
!===============

   Use OpenInit
   Use MLSCommon
   USE MLSCF
   USE MLSPCF3
   USE L3CF
   USE L3DZData
   USE L2GPData
   USE L2Interface
   USE OutputClose
   Use global_data
   Use DailyMapModule
   USE L3MMData
   USE L3MZData

   USE NumRecipesModule

   Implicit none

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- MonthlyCoreProcessing
!                
! This subroutine will process one product for all the input days data
!
! Remarks:  This module contains subroutines related to MLS L3 Monthly Map, Zonal Mean and Daily Zonal Mean Processing 
!

CONTAINS

!-------------------------------------------------------------------------
   SUBROUTINE MonthlyCoreProcessing(cfProd, pcf, cfDef, l2Days, l2gp, l3mm, mmA, mmD, &
					mzA, mzD, dzA, dzD, mis_l2Days, mis_Days)
!-------------------------------------------------------------------------

! Brief description of program
! This is the main program to run the Core processing.

! Parameters

	! Variable definitions

        TYPE( PCFMData_T ) :: pcf
        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T ) :: cfDef
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
	TYPE( L3DZData_T ), POINTER :: l3dz(:), dzA(:), dzD(:)

        TYPE( L3MMData_T ) :: l3mm, mmA, mmD
	TYPE( L3MZData_T ) :: l3mz, mzA, mzD

	Real, POINTER, DIMENSION(:, :, :) ::  alons(:, :, :), alats(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  dlons(:, :, :), dlats(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  atimes(:, :, :), dtimes(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  asolarTime(:, :, :), dsolarTime(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  asolarZenith(:, :, :), dsolarZenith(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)

	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 

        CHARACTER (LEN=480) :: msr

	INTEGER ::  error, l2Days, nlev, nf, nwv, numDays, pEndIndex, pStartIndex

        integer mis_l2Days

        CHARACTER (LEN=DATE_LEN) :: mis_Days(maxWindow)

        integer i, j, iP, kP, iD, iL
	real zonAvg, tau0, l3ret
	real lonD0_in, lonA0_in

!*** Initilize variables
 
	nlev = 0
        DO j = 1, l2gp(1)%nLevels
           IF( l2gp(1)%pressures(j) >= cfProd%l3presLvl(1) .AND. &
	       l2gp(1)%pressures(j) <= cfProd%l3presLvl(2) ) THEN
	     nlev = nlev + 1
	     IF (nlev == 1) pStartIndex = j 
	     pEndIndex = j 
	   ENDIF
        ENDDO

!*** Initilize POINTERS

        ALLOCATE( l3dz(l2Days), dzA(l2Days), dzD(l2Days), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3dz, l3r arrays.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

!!      Initialize Daily Zonal Mean 

        l3dz%name = cfProd%l3prodName
        dzA%name = TRIM(cfProd%l3prodName) // 'Ascending'
        dzD%name = TRIM(cfProd%l3prodName) // 'Descending'

        DO j = 1, l2Days
	  l3dz(j)%date = pcf%dates(j) 
	  dzA(j)%date = pcf%dates(j) 
	  dzD(j)%date = pcf%dates(j) 
        ENDDO

        DO j = 1, l2Days

           CALL AllocateL3DZ( nlev, cfDef%nNom, l3dz(j) )
           CALL AllocateL3DZ( nlev, cfDef%nNom, dzA(j) )
           CALL AllocateL3DZ( nlev, cfDef%nNom, dzD(j) )

           l3dz(j)%pressure = 0.0
           l3dz(j)%latitude = 0.0
           l3dz(j)%l3dzValue = 0.0
           l3dz(j)%l3dzPrecision = 0.0
           l3dz(j)%latRss = 0.0
           l3dz(j)%perMisPoints = 0
           dzA(j)%pressure = 0.0
           dzA(j)%latitude = 0.0
           dzA(j)%l3dzValue = 0.0
           dzA(j)%l3dzPrecision = 0.0
           dzA(j)%latRss = 0.0
           dzA(j)%perMisPoints = 0
           dzD(j)%pressure = 0.0
           dzD(j)%latitude = 0.0
           dzD(j)%l3dzValue = 0.0
           dzD(j)%l3dzPrecision = 0.0
           dzD(j)%latRss = 0.0
           dzD(j)%perMisPoints = 0

        ENDDO

!!      Initialize Monthly Zonal Mean 

        l3mz%name = cfProd%l3prodName
        mzA%name = TRIM(cfProd%l3prodName) // 'Ascending'
        mzD%name = TRIM(cfProd%l3prodName) // 'Descending'

        l3mz%startTime = l2gp(1)%time(1)
        l3mz%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)
        mzA%startTime = l2gp(1)%time(1)
        mzA%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)
        mzD%startTime = l2gp(1)%time(1)
        mzD%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)

        CALL AllocateL3MZ( nlev, cfDef%nNom, l3mz )
        CALL AllocateL3MZ( nlev, cfDef%nNom, mzA )
        CALL AllocateL3MZ( nlev, cfDef%nNom, mzD )

        l3mz%pressure = 0.0
        l3mz%latitude = 0.0
        l3mz%l3mzValue = 0.0
        l3mz%l3mzPrecision = 0.0
        l3mz%latRss = 0.0
        l3mz%perMisPoints = 0

        mzA%pressure = 0.0
        mzA%latitude = 0.0
        mzA%l3mzValue = 0.0
        mzA%l3mzPrecision = 0.0
        mzA%latRss = 0.0
        mzA%perMisPoints = 0

        mzD%pressure = 0.0
        mzD%latitude = 0.0
        mzD%l3mzValue = 0.0
        mzD%l3mzPrecision = 0.0
        mzD%latRss = 0.0
        mzD%perMisPoints = 0

!!      Initialize Monthly Map 

        l3mm%name = cfProd%l3prodName
        mmA%name = TRIM(cfProd%l3prodName) // 'Ascending'
        mmD%name = TRIM(cfProd%l3prodName) // 'Descending'

        l3mm%startTime = l2gp(1)%time(1)
        l3mm%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)
        mmA%startTime = l2gp(1)%time(1)
        mmA%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)
        mmD%startTime = l2gp(1)%time(1)
        mmD%endTime = l2gp(l2Days)%time(l2gp(l2Days)%nTimes)

        CALL AllocateL3MM( nlev, cfProd%nLats, cfProd%nLons, l3mm )
        CALL AllocateL3MM( nlev, cfProd%nLats, cfProd%nLons, mmA )
        CALL AllocateL3MM( nlev, cfProd%nLats, cfProd%nLons, mmD )

        l3mm%pressure = 0.0
        l3mm%latitude = 0.0
        l3mm%longitude = 0.0
        l3mm%l3mmValue = 0.0
        l3mm%l3mmPrecision = 0.0
        l3mm%perMisPoints = 0
        l3mm%misDays = mis_Days 

        mmA%pressure = 0.0
        mmA%latitude = 0.0
        mmA%longitude = 0.0
        mmA%l3mmValue = 0.0
        mmA%l3mmPrecision = 0.0
        mmA%perMisPoints = 0
        mmA%misDays = mis_Days 

        mmD%pressure = 0.0
        mmD%latitude = 0.0
        mmD%longitude = 0.0
        mmD%l3mmValue = 0.0
        mmD%l3mmPrecision = 0.0
        mmD%perMisPoints = 0
        mmD%misDays = mis_Days 

!*** Sort & Prepare the Data 

	Call SortDataMonthly(cfProd, l2Days, l2gp, 		&
			pStartIndex, pEndIndex,			&
			tau0, 					&
			anlats, dnlats, 			&
		        alats, dlats, 				&
			alons, dlons, 				&
			atimes, dtimes, 			&
			afields, dfields,			&
			asolarTime, asolarZenith,		&
			dsolarTime, dsolarZenith )

!*** Calculate Daily Zonal Means  

        Call DailyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days,        	&
                                    pStartIndex, pEndIndex,     &
                                    l3dz, dzA, dzD )

!*** Calculate Monthly Zonal Means  

        Call MonthlyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days,      	&
                                    pStartIndex, pEndIndex,     &
                                    l3mz, mzA, mzD )

!*** Calculate Monthly Map 

   	!Call MonthlyMapFromL2(cfProd, l2gp, 		&
	!		    pStartIndex, pEndIndex,  	&
	!		    anlats, dnlats, 		&
	!	            alats, dlats, 		&
	!		    alons, dlons, 		&
	!		    atimes, dtimes, 		&
	!		    afields, dfields, 		&
	!		    l3mm, mmA, mmD )
   	Call MonthlyMapFromL2_Simple(cfProd, cfDef, l2gp, l2Days, &
			    pStartIndex, pEndIndex, l3mm, mmA, mmD )

!*** Calculate Monthly Solar Parameter Ranges  

!   	Call CalSolarRange(cfProd, l2gp, 		&
!			    pStartIndex, pEndIndex,  	&
!			    anlats, dnlats, 		&
!		            asolarTime, dsolarTime,	&
!		            asolarZenith, dsolarZenith,	&
!			    l3mz, mzA, mzD )

!*** Allocate space for all the arrays

	DeAllocate(alats, dlats, alons, dlons, atimes, dtimes, afields, dfields)
	DeAllocate(asolarTime, asolarZenith, dsolarTime, dsolarZenith )

        CALL DestroyL3DZDatabase(l3dz)
        CALL DeallocateL3MZ(l3mz)

!-----------------------------------
   END SUBROUTINE MonthlyCoreProcessing
!-----------------------------------



!-------------------------------------------------------------------------
   SUBROUTINE SortDataMonthly(cfProd, l2Days, l2gp, pStartIndex, pEndIndex, tau0, anlats, dnlats, 	&
		alats_interp, dlats_interp, alons_interp, dlons_interp, 	&
		atimes_interp, dtimes_interp, afields_interp, dfields_interp, &
		asolarTime_interp, asolarZenith_interp, &
		dsolarTime_interp, dsolarZenith_interp)
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex, dindex_prev

        TYPE( L3CFMProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	Real, POINTER, DIMENSION(:, :) ::  l2Times(:, :), l2Lons_new(:, :), l2Lons(:, :),  l2Lons_old(:, :), &
					l2Lats(:, :), l2Values(:, :), l2SolarTime(:, :), l2SolarZenith(:, :)

	Real, POINTER, DIMENSION(:, :, :) ::  alons_interp(:, :, :), alats_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  dlons_interp(:, :, :), dlats_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  atimes_interp(:, :, :), dtimes_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  asolarTime_interp(:, :, :), dsolarTime_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  asolarZenith_interp(:, :, :), dsolarZenith_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  afields_interp(:, :, :), dfields_interp(:, :, :)

	Real dlons, alons, lons_found, lats_found, times_found, fields_found, 	&
		    solarTime_found, solarZenith_found, &
		  slope, slope_time, slope_solarTime, slope_solarZenith, slope_field, sTime, ddtad

	Real lons_found_old

	Real tau0
 
	Integer nPd, pStartIndex, pEndIndex, iMin

	INTEGER, DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: anlats, dnlats

!*** Calculate the number of points in each pressure level 

	nPd = 15

        nterms = 0
	DO I = 1, l2Days
            nterms = nterms + l2gp(I)%nTimes
	ENDDO

!*** Allocate space for all the arrays

        ALLOCATE(l2Times(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Lons(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Lats(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Values(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2SolarTime(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2SolarZenith(nterms, pEndIndex-pStartIndex+1), STAT=error)

        ALLOCATE(l2Lons_new(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Lons_old(nterms, pEndIndex-pStartIndex+1), STAT=error)

        ALLOCATE(alons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(alats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(atimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(asolarTime_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(asolarZenith_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(afields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dtimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dsolarTime_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dsolarZenith_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dfields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

!*** Re-arrange the data into longitude order for each pressure level 

        nstart = 1
	DO iD = 1, l2Days
	  DO iT = 1, l2gp(iD)%nTimes
	    iP = 0
	    DO kP = pStartIndex, pEndIndex 
	        iP = iP + 1
                l2Times(nstart, iP) = l2gp(iD)%time(iT)
                l2SolarTime(nstart, iP) = l2gp(iD)%solarTime(iT)
                l2SolarZenith(nstart, iP) = l2gp(iD)%solarZenith(iT)
                l2Lons(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lons_old(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lats(nstart, iP) = l2gp(iD)%latitude(iT)
                l2Values(nstart, iP) = l2gp(iD)%l2gpValue(1, kP, iT) 
	    ENDDO
            nstart = nstart + 1
	  ENDDO
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
		slope = (l2Lats(lindex+1, iP)-l2Lats(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
		slope_time = (l2Times(lindex+1, iP)-l2Times(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
		slope_field = (l2Values(lindex+1, iP)-l2Values(lindex, iP))/	&
			(l2Lons(lindex+1, iP)-l2Lons(lindex, iP))
		slope_solarTime = (l2SolarTime(lindex+1, iP)-l2SolarTime(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
		slope_solarZenith = (l2SolarZenith(lindex+1, iP)-l2SolarZenith(lindex, iP))/	&
			(l2Lons_new(lindex+1, iP)-l2Lons_new(lindex, iP))
		IF(ABS(slope) < 1.e-6) THEN
		   lons_found = l2Lons_new(lindex, iP)
		   lons_found_old = l2Lons_old(lindex, iP)
		ELSE
		   lons_found = l2Lons_new(lindex, iP) + (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		   lons_found_old = l2Lons_old(lindex, iP) + (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		END IF
		lats_found = cfProd%latGridMap(J)
		times_found = l2Times(lindex, iP) + (lons_found-l2Lons_new(lindex, iP))*slope_time
		solarTime_found = l2SolarTime(lindex, iP) + (lons_found-l2Lons_new(lindex, iP))*slope_solarTime
		solarZenith_found = l2SolarZenith(lindex, iP) + (lons_found-l2Lons_new(lindex, iP))*slope_solarZenith
		fields_found = l2Values(lindex, iP) + (lons_found-l2Lons_new(lindex, iP))*slope_field
                aindex = J
		EXIT
	      END IF
	     ENDDO
              
	     IF( slope < 0 .and. aindex /= aindex_prev) THEN
	         aindex_prev = aindex
		 !IF(dnlats(aindex, iP) > 0) THEN
       	           anlats(aindex, iP) = anlats(aindex, iP) + 1 
	           alons_interp(aindex, anlats(aindex, iP), iP) = lons_found
	           alats_interp(aindex, anlats(aindex, iP), iP) = lats_found
	           atimes_interp(aindex, anlats(aindex, iP), iP) = times_found
	           asolarTime_interp(aindex, anlats(aindex, iP), iP) = solarTime_found
	           asolarZenith_interp(aindex, anlats(aindex, iP), iP) = solarZenith_found
	           afields_interp(aindex, anlats(aindex, iP), iP) = fields_found
	         !END IF
	     END IF

	     IF( slope >= 0 .and. aindex /= dindex_prev) THEN
	         dindex_prev = aindex
		 IF(anlats(aindex, iP) > 0) THEN
       	           dnlats(aindex, iP) = dnlats(aindex, iP) + 1 
	           dlons_interp(aindex, dnlats(aindex, iP), iP) = lons_found
	           dlats_interp(aindex, dnlats(aindex, iP), iP) = lats_found
	           dtimes_interp(aindex, dnlats(aindex, iP), iP) = times_found
	           dsolarTime_interp(aindex, dnlats(aindex, iP), iP) = solarTime_found
	           dsolarZenith_interp(aindex, dnlats(aindex, iP), iP) = solarZenith_found
	           dfields_interp(aindex, dnlats(aindex, iP), iP) = fields_found
	         END IF
	     END IF


	    END IF

          ENDDO


        ENDDO

!*** Deallocate intermidiate arrays 

        DEALLOCATE(l2Times, l2SolarTime, l2SolarZenith, l2Lons, l2Lons_new, l2Lons_old, l2Lats, l2Values)

!-----------------------------------
   END SUBROUTINE SortDataMonthly
!-----------------------------------



!-------------------------------------------------------------------------
   SUBROUTINE DailyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days, pStartIndex, pEndIndex,  	&
			     l3dz, dzA, dzD )
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr
	INTEGER ::  iCom, iAasc, iDes, nloop, aindex

        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T ) :: cfDef
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
	TYPE( L3DZData_T ), POINTER :: l3dz(:), dzA(:), dzD(:)

	Integer, POINTER, DIMENSION(:) ::  iComArr(:), iAscArr(:), iDesArr(:)
	Real, POINTER, DIMENSION(:,:) ::   comFieldArr(:, :), ascFieldArr(:, :), desFieldArr(:, :)

	Integer l2Days, pStartIndex, pEndIndex, iMin

	Real slope

!*** Allocate space for all the arrays

        ALLOCATE(iComArr(cfDef%nNom), STAT=error)
        ALLOCATE(iAscArr(cfDef%nNom), STAT=error)
        ALLOCATE(iDesArr(cfDef%nNom), STAT=error)

        ALLOCATE(comFieldArr(cfDef%nNom, 100), STAT=error)
        ALLOCATE(ascFieldArr(cfDef%nNom, 100), STAT=error)
        ALLOCATE(desFieldArr(cfDef%nNom, 100), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

	DO iT = 1, cfDef%nNom 
	  iComArr(iT) = 0
	  iAscArr(iT) = 0
	  iDesArr(iT) = 0
	  DO J = 1, 100 
	  	comFieldArr(iT, J) = 0.0
	  	ascFieldArr(iT, J) = 0.0
	  	desFieldArr(iT, J) = 0.0
          ENDDO
        ENDDO

!*** initialization 

        DO iD = 1, l2Days
	  iP = 0
	  DO kP = pStartIndex, pEndIndex 
	    iP = iP + 1
	    l3dz(iD)%pressure(iP) = l2gp(1)%pressures(kP) 
	    dzA(iD)%pressure(iP)  = l2gp(1)%pressures(kP) 
	    dzD(iD)%pressure(iP)  = l2gp(1)%pressures(kP) 
          ENDDO

          DO J = 1, cfDef%nNom 
	    l3dz(iD)%latitude(J) = cfDef%l2nomLats(J) 
	    dzA(iD)%latitude(J)  = cfDef%l2nomLats(J) 
	    dzD(iD)%latitude(J)  = cfDef%l2nomLats(J) 
	  END DO
        ENDDO

        DO iD = 1, l2Days
          DO J = 1, cfDef%nNom 
	    l3dz(iD)%localSolarZenithAngle(1, J) =  1.e20 
	    l3dz(iD)%localSolarZenithAngle(2, J) = -1.e20 
	    l3dz(iD)%localSolarTime(1, J) =  1.e20 
	    l3dz(iD)%localSolarTime(2, J) = -1.e20 
            DO I = 1, l3dz(iD)%nLevels 
	    	l3dz(iD)%perMisPoints(I,J) = 0 
            ENDDO

	    dzA(iD)%localSolarZenithAngle(1, J) =  1.e20 
	    dzA(iD)%localSolarZenithAngle(2, J) = -1.e20 
	    dzA(iD)%localSolarTime(1, J) =  1.e20 
	    dzA(iD)%localSolarTime(2, J) = -1.e20 
            DO I = 1, dzA(iD)%nLevels 
	    	dzA(iD)%perMisPoints(I,J) = 0 
            ENDDO

	    dzD(iD)%localSolarZenithAngle(1, J) =  1.e20 
	    dzD(iD)%localSolarZenithAngle(2, J) = -1.e20 
	    dzD(iD)%localSolarTime(1, J) =  1.e20 
	    dzD(iD)%localSolarTime(2, J) = -1.e20 
            DO I = 1, dzD(iD)%nLevels 
	    	dzD(iD)%perMisPoints(I,J) = 0 
            ENDDO
          ENDDO
        ENDDO

!*** Re-arrange the data into longitude order for each pressure level 

        DO iD = 1, l2Days
	  iP = 0
	  DO kP = pStartIndex, pEndIndex 
	    iP = iP + 1

	!** Initilizing
	    DO iT = 1, cfDef%nNom 
	      iComArr(iT) = 0
	      iAscArr(iT) = 0
	      iDesArr(iT) = 0

	      l3dz(iD)%l3dzValue(iP, iT) = 0.0 
	      dzA(iD)%l3dzValue(iP, iT) = 0.0 
	      dzD(iD)%l3dzValue(iP, iT) = 0.0 
            ENDDO

	    DO iT = 1, l2gp(iD)%nTimes-1
   	      iCom = FindIndexForNormGrid(cfDef, l2gp(iD)%latitude(iT))
	      l3dz(iD)%l3dzValue(iP, iCom) = l3dz(iD)%l3dzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
      !** Solar Zenith Angle & Time
	      if(iP == 1) then
	         if(l3dz(iD)%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	            l3dz(iD)%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		 else if(l3dz(iD)%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	            l3dz(iD)%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		 end if

	         if(l3dz(iD)%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	            l3dz(iD)%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		 else if(l3dz(iD)%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	            l3dz(iD)%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		 end if
	      end if
	      iComArr(iCom) = iComArr(iCom) + 1 
	      comFieldArr(iCom, iComArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	      IF( l2gp(iD)%longitude(iT) >= -PI .AND. &
		  l2gp(iD)%longitude(iT) < 0.0  .AND. &
		  l2gp(iD)%longitude(iT+1) <= PI .AND. &
		  l2gp(iD)%longitude(iT+1) > 0.0 .AND. &
		  l2gp(iD)%longitude(iT+1) > l2gp(iD)%longitude(iT) ) THEN
		slope = (l2gp(iD)%longitude(iT+1)-2.0*PI-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	      ELSE
		slope = (l2gp(iD)%longitude(iT+1)-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	      END IF
	      IF( slope < 0.0 ) THEN
	        dzA(iD)%l3dzValue(iP, iCom) = dzA(iD)%l3dzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iAscArr(iCom) = iAscArr(iCom) + 1 
	        ascFieldArr(iCom, iAscArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        if(iP == 1) then
	           if(dzA(iD)%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	              dzA(iD)%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		   else if(dzA(iD)%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	              dzA(iD)%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		   end if

	           if(dzA(iD)%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	              dzA(iD)%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		   else if(dzA(iD)%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	              dzA(iD)%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		   end if
	        end if
	      ELSE
	        dzD(iD)%l3dzValue(iP, iCom) = dzD(iD)%l3dzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iDesArr(iCom) = iDesArr(iCom) + 1 
	        desFieldArr(iCom, iDesArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        if(iP == 1) then
	           if(dzD(iD)%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	              dzD(iD)%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		   else if(dzD(iD)%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	              dzD(iD)%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		   end if

	           if(dzD(iD)%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	              dzD(iD)%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		   else if(dzD(iD)%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	              dzD(iD)%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		   end if
	        end if
	      END IF	
	    ENDDO

	    DO iT = 1, cfDef%nNom 
	      IF(iComArr(iT) == 0) THEN
	      	l3dz(iD)%l3dzValue(iP, iT) = 0.0 
	      	dzA(iD)%l3dzValue(iP, iT) = 0.0 
	      	dzD(iD)%l3dzValue(iP, iT) = 0.0 

	      	l3dz(iD)%latRss(iP, iT) = 0.0 
	      	dzA(iD)%latRss(iP, iT) = 0.0 
	      	dzD(iD)%latRss(iP, iT) = 0.0 
	      ELSE
	      	l3dz(iD)%l3dzValue(iP, iT) = l3dz(iD)%l3dzValue(iP, iT)/real(iComArr(iT))
		IF(iAscArr(iT) > 0) THEN
	      		dzA(iD)%l3dzValue(iP, iT) = dzA(iD)%l3dzValue(iP, iT)/real(iAscArr(iT))
		ELSE
	      		dzA(iD)%l3dzValue(iP, iT) = 0.0 
		END IF
		IF(iDesArr(iT) > 0) THEN
	      		dzD(iD)%l3dzValue(iP, iT) = dzD(iD)%l3dzValue(iP, iT)/real(iDesArr(iT))
		ELSE
	      		dzD(iD)%l3dzValue(iP, iT) = 0.0 
		END IF

		!*** calculate Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

	  	DO J = 1, iComArr(iT)  
	  		l3dz(iD)%latRss(iP, iT) = l3dz(iD)%latRss(iP, iT) + 				&
						  (comFieldArr(iT, J)-l3dz(iD)%l3dzValue(iP, iT))* 	&
						  (comFieldArr(iT, J)-l3dz(iD)%l3dzValue(iP, iT))
          	ENDDO
	  	l3dz(iD)%latRss(iP, iT) = l3dz(iD)%latRss(iP, iT)/real(iComArr(iT))

		IF(iAscArr(iT) > 0) THEN
	  	   DO J = 1, iAscArr(iT)  
	  		dzA(iD)%latRss(iP, iT) = dzA(iD)%latRss(iP, iT) + 				&
						 (ascFieldArr(iT, J)-dzA(iD)%l3dzValue(iP, iT))* 	&
						 (ascFieldArr(iT, J)-dzA(iD)%l3dzValue(iP, iT))
          	   ENDDO
	  	   dzA(iD)%latRss(iP, iT) = dzA(iD)%latRss(iP, iT)/real(iAscArr(iT))
		ELSE
	  	   dzA(iD)%latRss(iP, iT) = 0.0 
		END IF

		IF(iDesArr(iT) > 0) THEN
	  	   DO J = 1, iDesArr(iT)  
	  		dzD(iD)%latRss(iP, iT) = dzD(iD)%latRss(iP, iT) + 				&
						 (desFieldArr(iT, J)-dzD(iD)%l3dzValue(iP, iT))* 	&
						 (desFieldArr(iT, J)-dzD(iD)%l3dzValue(iP, iT))
          	   ENDDO
	  	   dzD(iD)%latRss(iP, iT) = dzD(iD)%latRss(iP, iT)/real(iDesArr(iT))
		ELSE
	  	   dzD(iD)%latRss(iP, iT) = 0.0 
		END IF

	      END IF
	    ENDDO

	  ENDDO
	ENDDO

!*** Deallocate 

        DEALLOCATE(iComArr, iAscArr, iDesArr)
        DEALLOCATE(comFieldArr, ascFieldArr, desFieldArr)

!-----------------------------------
   END SUBROUTINE DailyZonalMeanFromL2
!-----------------------------------

!-------------------------------------------------------------------------
   SUBROUTINE MonthlyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days, pStartIndex, pEndIndex,  	&
			     	l3mz, mzA, mzD )
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr
	INTEGER ::  iCom, iAsc, iDes, nloop, aindex

        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T ) :: cfDef
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
	TYPE( L3MZData_T ) :: l3mz, mzA, mzD

	Integer, POINTER, DIMENSION(:) ::  iComArr(:), iAscArr(:), iDesArr(:)
	Real, POINTER, DIMENSION(:,:) ::   comFieldArr(:, :), ascFieldArr(:, :), desFieldArr(:, :)

	Integer l2Days, pStartIndex, pEndIndex, iMin

	Real slope

!*** Allocate space for all the arrays

        ALLOCATE(iComArr(cfDef%nNom), STAT=error)
        ALLOCATE(iAscArr(cfDef%nNom), STAT=error)
        ALLOCATE(iDesArr(cfDef%nNom), STAT=error)

        ALLOCATE(comFieldArr(cfDef%nNom, 100*l2Days), STAT=error)
        ALLOCATE(ascFieldArr(cfDef%nNom, 100*l2Days), STAT=error)
        ALLOCATE(desFieldArr(cfDef%nNom, 100*l2Days), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

!*** initialization 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  l3mz%pressure(iP) = l2gp(1)%pressures(kP) 
	  mzA%pressure(iP)  = l2gp(1)%pressures(kP) 
	  mzD%pressure(iP)  = l2gp(1)%pressures(kP) 
          DO J = 1, cfDef%nNom 
	      l3mz%l3mzValue(iP, J) = 0.0 
	      mzA%l3mzValue(iP, J) = 0.0 
	      mzD%l3mzValue(iP, J) = 0.0 
          ENDDO
        ENDDO

        DO J = 1, cfDef%nNom 
	     l3mz%latitude(J) = cfDef%l2nomLats(J) 
	     mzA%latitude(J)  = cfDef%l2nomLats(J) 
	     mzD%latitude(J)  = cfDef%l2nomLats(J) 
	END DO

        DO J = 1, cfDef%nNom 
	    l3mz%localSolarZenithAngle(1, J) =  1.e20 
	    l3mz%localSolarZenithAngle(2, J) = -1.e20 
	    l3mz%localSolarTime(1, J) =  1.e20 
	    l3mz%localSolarTime(2, J) = -1.e20 
            DO I = 1, l3mz%nLevels 
	    	l3mz%perMisPoints(I,J) = 0 
            ENDDO

	    mzA%localSolarZenithAngle(1, J) =  1.e20 
	    mzA%localSolarZenithAngle(2, J) = -1.e20 
	    mzA%localSolarTime(1, J) =  1.e20 
	    mzA%localSolarTime(2, J) = -1.e20 
            DO I = 1, mzA%nLevels 
	    	mzA%perMisPoints(I,J) = 0 
            ENDDO

	    mzD%localSolarZenithAngle(1, J) =  1.e20 
	    mzD%localSolarZenithAngle(2, J) = -1.e20 
	    mzD%localSolarTime(1, J) =  1.e20 
	    mzD%localSolarTime(2, J) = -1.e20 
            DO I = 1, mzD%nLevels 
	    	mzD%perMisPoints(I,J) = 0 
            ENDDO
        ENDDO

!*** Re-arrange the data into longitude order for each pressure level 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1

	  DO iT = 1, cfDef%nNom 
	    iComArr(iT) = 0
	    iAscArr(iT) = 0
	    iDesArr(iT) = 0
	    DO J = 1,  100*l2Days
	  	comFieldArr(iT, J) = 0.0
	  	ascFieldArr(iT, J) = 0.0
	  	desFieldArr(iT, J) = 0.0
            ENDDO
          ENDDO

          DO iD = 1, l2Days
	    DO iT = 1, l2gp(iD)%nTimes-1
	      !iCom = real(l2gp(iD)%latitude(iT)-l3mz%latitude(1))/real(l3mz%latitude(2)-l3mz%latitude(1))+1.5
   	      iCom = FindIndexForNormGrid(cfDef, l2gp(iD)%latitude(iT))
	      l3mz%l3mzValue(iP, iCom) = l3mz%l3mzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
      !** Solar Zenith Angle & Time
	      if(iP == 1) then
	         if(l3mz%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	            l3mz%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		 else if(l3mz%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	            l3mz%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		 end if

	         if(l3mz%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	            l3mz%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		 else if(l3mz%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	            l3mz%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		 end if
	      end if
	      iComArr(iCom) = iComArr(iCom) + 1 
	      comFieldArr(iCom, iComArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	      IF( l2gp(iD)%longitude(iT) >= -PI .AND. &
		  l2gp(iD)%longitude(iT) < 0.0  .AND. &
		  l2gp(iD)%longitude(iT+1) <= PI .AND. &
		  l2gp(iD)%longitude(iT+1) > 0.0 .AND. &
		  l2gp(iD)%longitude(iT+1) > l2gp(iD)%longitude(iT) ) THEN
		slope = (l2gp(iD)%longitude(iT+1)-2.0*PI-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	      ELSE
		slope = (l2gp(iD)%longitude(iT+1)-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	      END IF
	      IF( slope < 0.0 ) THEN
	        mzA%l3mzValue(iP, iCom) = mzA%l3mzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iAscArr(iCom) = iAscArr(iCom) + 1 
	        ascFieldArr(iCom, iAscArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        if(iP == 1) then
	           if(mzA%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	              mzA%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		   else if(mzA%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	              mzA%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		   end if

	           if(mzA%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	              mzA%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		   else if(mzA%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	              mzA%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		   end if
	        end if
	      ELSE
	        mzD%l3mzValue(iP, iCom) = mzD%l3mzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iDesArr(iCom) = iDesArr(iCom) + 1 
	        desFieldArr(iCom, iDesArr(iCom)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        if(iP == 1) then
	           if(mzD%localSolarZenithAngle(2, iCom) <= l2gp(iD)%solarZenith(iT)) then
	              mzD%localSolarZenithAngle(2, iCom) =  l2gp(iD)%solarZenith(iT)
		   else if(mzD%localSolarZenithAngle(1, iCom) >= l2gp(iD)%solarZenith(iT)) then
	              mzD%localSolarZenithAngle(1, iCom) =  l2gp(iD)%solarZenith(iT)
		   end if

	           if(mzD%localSolarTime(2, iCom) <= l2gp(iD)%solarTime(iT)) then
	              mzD%localSolarTime(2, iCom) =  l2gp(iD)%solarTime(iT)
		   else if(mzD%localSolarTime(1, iCom) >= l2gp(iD)%solarTime(iT)) then
	              mzD%localSolarTime(1, iCom) =  l2gp(iD)%solarTime(iT)
		   end if
	        end if
	      END IF	
	    ENDDO
	  ENDDO

	  DO iT = 1, cfDef%nNom 
	    IF(iComArr(iT) == 0) THEN
	      	l3mz%l3mzValue(iP, iT) = 0.0 
	      	mzA%l3mzValue(iP, iT) = 0.0 
	      	mzD%l3mzValue(iP, iT) = 0.0 

	      	l3mz%latRss(iP, iT) = 0.0 
	      	mzA%latRss(iP, iT) = 0.0 
	      	mzD%latRss(iP, iT) = 0.0 
	    ELSE
	        l3mz%l3mzValue(iP, iT) = l3mz%l3mzValue(iP, iT)/real(iComArr(iT))
		IF(iAscArr(iT) > 0) THEN
	    		mzA%l3mzValue(iP, iT) = mzA%l3mzValue(iP, iT)/real(iAscArr(iT))
		ELSE
	    		mzA%l3mzValue(iP, iT) = 0.0 
		END IF
		IF(iDesArr(iT) > 0) THEN
	    		mzD%l3mzValue(iP, iT) = mzD%l3mzValue(iP, iT)/real(iDesArr(iT))
		ELSE
	    		mzD%l3mzValue(iP, iT) = 0.0 
		END IF

		!*** calculate Root-Sum-Square for each latitude, dimensioned (nLevels, nLats)

	  	DO J = 1, iComArr(iT)  
	  		l3mz%latRss(iP, iT) = l3mz%latRss(iP, iT) + 				&
					    	  (comFieldArr(iT, J)-l3mz%l3mzValue(iP, iT))* 	&
						  (comFieldArr(iT, J)-l3mz%l3mzValue(iP, iT))
          	ENDDO
	  	l3mz%latRss(iP, iT) = l3mz%latRss(iP, iT)/real(iComArr(iT))

		IF(iAscArr(iT) > 0) THEN
	  	   DO J = 1, iAscArr(iT)  
	  		mzA%latRss(iP, iT) = mzA%latRss(iP, iT) + 				&
						 (ascFieldArr(iT, J)-mzA%l3mzValue(iP, iT))* 	&
						 (ascFieldArr(iT, J)-mzA%l3mzValue(iP, iT))
          	   ENDDO
	  	   mzA%latRss(iP, iT) = mzA%latRss(iP, iT)/real(iAscArr(iT))
		ELSE
	  	   mzA%latRss(iP, iT) = 0.0 
		END IF

		IF(iDesArr(iT) > 0) THEN
	  	   DO J = 1, iDesArr(iT)  
	  		mzD%latRss(iP, iT) = mzD%latRss(iP, iT) + 				&
						 (desFieldArr(iT, J)-mzD%l3mzValue(iP, iT))* 	&
						 (desFieldArr(iT, J)-mzD%l3mzValue(iP, iT))
          	   ENDDO
	  	   mzD%latRss(iP, iT) = mzD%latRss(iP, iT)/real(iDesArr(iT))
		ELSE
	  	   mzD%latRss(iP, iT) = 0.0 
		END IF

	    END IF
	  ENDDO

	ENDDO

!*** Deallocate 

        DEALLOCATE(iComArr, iAscArr, iDesArr)
        DEALLOCATE(comFieldArr, ascFieldArr, desFieldArr)

!-----------------------------------
   END SUBROUTINE MonthlyZonalMeanFromL2
!-----------------------------------

!-------------------------------------------------------------------------
   SUBROUTINE MonthlyMapFromL2_Simple(cfProd, cfDef, l2gp, l2Days, pStartIndex, pEndIndex,  	&
			     	l3mm, mmA, mmD )
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, kP, nterms, nstart, nr
	INTEGER ::  iComLat, iComLon, iAsc, iDes, nloop, aindex

        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T )  :: cfDef
        TYPE( L3MMData_T )  :: l3mm, mmA, mmD

        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	Integer, POINTER, DIMENSION(:,:) ::  iComArr(:,:), iAscArr(:,:), iDesArr(:,:)
	Real, POINTER, DIMENSION(:,:,:) ::   comFieldArr(:,:,:), ascFieldArr(:,:,:), desFieldArr(:,:,:)

	Integer l2Days, pStartIndex, pEndIndex, nc

	Real slope, lonStart, newField

!*** Allocate space for all the arrays

        ALLOCATE(iComArr(cfProd%nLats,cfProd%nLons), STAT=error)
        ALLOCATE(iAscArr(cfProd%nLats,cfProd%nLons), STAT=error)
        ALLOCATE(iDesArr(cfProd%nLats,cfProd%nLons), STAT=error)

        ALLOCATE(comFieldArr(cfProd%nLats,cfProd%nLons, 10*l2Days), STAT=error)
        ALLOCATE(ascFieldArr(cfProd%nLats,cfProd%nLons, 10*l2Days), STAT=error)
        ALLOCATE(desFieldArr(cfProd%nLats,cfProd%nLons, 10*l2Days), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

!*** initialization 

        iP = 0
	DO kP = pStartIndex, pEndIndex 
          iP = iP + 1
	  l3mm%pressure(iP) = l2gp(1)%pressures(kP) 
	  mmA%pressure(iP)  = l2gp(1)%pressures(kP) 
	  mmD%pressure(iP)  = l2gp(1)%pressures(kP) 
	  DO I = 1, cfProd%nLats 
	  DO J = 1, cfProd%nLons 
	      l3mm%l3mmValue(iP, I, J) = 0.0 
	      mmA%l3mmValue(iP, I, J) = 0.0 
	      mmD%l3mmValue(iP, I, J) = 0.0 
          ENDDO
          ENDDO
	END DO

        DO J = 1, cfProd%nLats
	     l3mm%latitude(J) = cfProd%latGridMap(J) 
	     mmA%latitude(J)  = cfProd%latGridMap(J) 
	     mmD%latitude(J)  = cfProd%latGridMap(J) 
	END DO

        DO K = 1, cfProd%nLons
	     l3mm%longitude(K) = cfProd%longGrid(K) 
	     mmA%longitude(K)  = cfProd%longGrid(K) 
	     mmD%longitude(K)  = cfProd%longGrid(K) 
	END DO

        DO I = 1, l3mm%nLevels 
	    	l3mm%perMisPoints(I) = 0 
        ENDDO
        DO I = 1, mmA%nLevels 
	    	mmA%perMisPoints(I) = 0 
        ENDDO
        DO I = 1, mmD%nLevels 
	    	mmD%perMisPoints(I) = 0 
        ENDDO

!*** Start calculation 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1

	  DO I = 1, cfProd%nLats 
	  DO J = 1, cfProd%nLons 
	    iComArr(I,J) = 0
	    iAscArr(I,J) = 0
	    iDesArr(I,J) = 0
	    DO K = 1,  10*l2Days
	  	comFieldArr(I,J, K) = 0.0
	  	ascFieldArr(I,J, K) = 0.0
	  	desFieldArr(I,J, K) = 0.0
            ENDDO
          ENDDO
          ENDDO

          DO iD = 1, l2Days
	    DO iT = 1, l2gp(iD)%nTimes-1
   	      iComLat = FindLatIndexForL3Grid(cfProd, l2gp(iD)%latitude(iT))
   	      iComLon = FindLonIndexForL3Grid(cfProd, l2gp(iD)%longitude(iT))
	      IF(l2gp(iD)%l2gpValue(1, kP, iT) > 0.0) THEN
	        l3mm%l3mmValue(iP, iComLat, iComLon) = l3mm%l3mmValue(iP, iComLat, iComLon) + l2gp(iD)%l2gpValue(1, kP, iT)
	        l3mm%l3mmPrecision(iP, iComLat, iComLon) = 0.0 
	        iComArr(iComLat, iComLon) = iComArr(iComLat, iComLon) + 1 
	        comFieldArr(iComLat, iComLon, iComArr(iComLat, iComLon)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        IF( l2gp(iD)%longitude(iT) >= -PI .AND. &
		  l2gp(iD)%longitude(iT) < 0.0  .AND. &
		  l2gp(iD)%longitude(iT+1) <= PI .AND. &
		  l2gp(iD)%longitude(iT+1) > 0.0 .AND. &
		  l2gp(iD)%longitude(iT+1) > l2gp(iD)%longitude(iT) ) THEN
		  slope = (l2gp(iD)%longitude(iT+1)-2.0*PI-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	        ELSE
		  slope = (l2gp(iD)%longitude(iT+1)-l2gp(iD)%longitude(iT))/	&
		        (l2gp(iD)%latitude(iT+1)-l2gp(iD)%latitude(iT)) 
 	        END IF
	        IF( slope < 0.0 ) THEN
	          mmA%l3mmValue(iP, iComLat, iComLon) = mmA%l3mmValue(iP, iComLat, iComLon) + l2gp(iD)%l2gpValue(1, kP, iT)
	          iAscArr(iComLat, iComLon) = iAscArr(iComLat, iComLon) + 1 
	          ascFieldArr(iComLat, iComLon, iAscArr(iComLat, iComLon)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        ELSE
	          mmD%l3mmValue(iP, iComLat, iComLon) = mmD%l3mmValue(iP, iComLat, iComLon) + l2gp(iD)%l2gpValue(1, kP, iT)
	          iDesArr(iComLat, iComLon) = iDesArr(iComLat, iComLon) + 1 
	          desFieldArr(iComLat, iComLon, iDesArr(iComLat, iComLon)) = l2gp(iD)%l2gpValue(1, kP, iT)
	        END IF	
	      END IF	
	    ENDDO
	  ENDDO

	  DO I = 1, cfProd%nLats 
	  DO J = 1, cfProd%nLons 
	    IF(iComArr(I, J) == 0) THEN
	      	l3mm%l3mmValue(iP, I, J) = 0.0 
	      	mmA%l3mmValue(iP, I, J) = 0.0 
	      	mmD%l3mmValue(iP, I, J) = 0.0 
	    ELSE
	        l3mm%l3mmValue(iP, I, J) = l3mm%l3mmValue(iP, I, J)/real(iComArr(I, J))
		IF(iAscArr(I, J) > 0) THEN
	    		mmA%l3mmValue(iP, I, J) = mmA%l3mmValue(iP, I, J)/real(iAscArr(I, J))
		ELSE
	    		mmA%l3mmValue(iP, I, J) = 0.0 
		END IF
		IF(iDesArr(I, J) > 0) THEN
	    		mmD%l3mmValue(iP, I, J) = mmD%l3mmValue(iP, I, J)/real(iDesArr(I, J))
		ELSE
	    		mmD%l3mmValue(iP, I, J) = 0.0 
		END IF
	    END IF
	  ENDDO
	  ENDDO

        ENDDO

!*** Deallocate 

        DEALLOCATE(iComArr, iAscArr, iDesArr)
        DEALLOCATE(comFieldArr, ascFieldArr, desFieldArr)

!-----------------------------------
   END SUBROUTINE MonthlyMapFromL2_Simple
!-----------------------------------


!-------------------------------------------------------------------------
   SUBROUTINE MonthlyMapFromL2(cfProd, l2gp, pStartIndex, pEndIndex,  	&
			    	anlats, dnlats, alats, dlats, alons, dlons, atimes, dtimes, afields, dfields, 	&
			     	l3mm, mmA, mmD )
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, kP, nterms, nstart, nr

        TYPE( L3CFMProd_T ) :: cfProd
        TYPE( L3MMData_T ) :: l3mm, mmA, mmD

        TYPE( L2GPData_T ), POINTER :: l2gp(:)
	Real, POINTER, DIMENSION(:, :, :) ::  alons(:, :, :), alats(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  dlons(:, :, :), dlats(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  atimes(:, :, :), dtimes(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)

	Real, POINTER, DIMENSION(:) ::  lons_rev(:), fields_rev(:), Y2(:)

	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 

	Integer pStartIndex, pEndIndex, nc

	Real slope, lonStart, newField

!*** initialization 

        iP = 0
	DO kP = pStartIndex, pEndIndex 
          iP = iP + 1
	  l3mm%pressure(iP) = l2gp(1)%pressures(kP) 
	  mmA%pressure(iP)  = l2gp(1)%pressures(kP) 
	  mmD%pressure(iP)  = l2gp(1)%pressures(kP) 
	  DO I = 1, cfProd%nLats 
	  DO J = 1, cfProd%nLons 
	      l3mm%l3mmValue(iP, I, J) = 0.0 
	      mmA%l3mmValue(iP, I, J) = 0.0 
	      mmD%l3mmValue(iP, I, J) = 0.0 
          ENDDO
          ENDDO
	END DO

        DO J = 1, cfProd%nLats
	     l3mm%latitude(J) = cfProd%latGridMap(J) 
	     mmA%latitude(J)  = cfProd%latGridMap(J) 
	     mmD%latitude(J)  = cfProd%latGridMap(J) 
	END DO

        DO K = 1, cfProd%nLons
	     l3mm%longitude(K) = cfProd%longGrid(K) 
	     mmA%longitude(K)  = cfProd%longGrid(K) 
	     mmD%longitude(K)  = cfProd%longGrid(K) 
	END DO

        DO I = 1, l3mm%nLevels 
	    	l3mm%perMisPoints(I) = 0 
        ENDDO
        DO I = 1, mmA%nLevels 
	    	mmA%perMisPoints(I) = 0 
        ENDDO
        DO I = 1, mmD%nLevels 
	    	mmD%perMisPoints(I) = 0 
        ENDDO

!*** Start calculation 

        iP = 0
	DO kP = pStartIndex, pEndIndex 
          iP = iP + 1

          DO J = 1, cfProd%nLats

	     ! Ascending Mode

             IF( anlats(J, iP) > 0 ) THEN

!*** Allocate space for all the arrays
        	ALLOCATE(lons_rev(anlats(J, iP)), STAT=error)
        	ALLOCATE(fields_rev(anlats(J, iP)), STAT=error)

                DO i = 1, anlats(J, iP)
                    alons(J, i, iP) = alons(J, i, iP) + (i-1)*2.0*PI
                    lons_rev(i) = alons(J, anlats(J, iP)+1-i, iP)
                    fields_rev(i) = afields(J, anlats(J, iP)+1-i, iP)
                ENDDO 
                ALLOCATE(Y2(anlats(J, iP)), STAT=error)
		CALL SPLINE(lons_rev, fields_rev, anlats(J, iP), 0.0, 0.0, Y2)
                DO K = 1, cfProd%nLons
		  IF( cfProd%longGrid(K) > lons_rev(1) ) THEN
		     lonStart = cfProd%longGrid(K) - 2.0*PI
		  ELSE
		     lonStart = cfProd%longGrid(K)
		  END IF
		  nc = 0
		  mmA%l3mmValue(iP, J, K) = 0.0 
		  DO
		     nc = nc + 1
		     CALL SPLINT(lons_rev, fields_rev, Y2, anlats(J, iP), lonStart, newField)
		     mmA%l3mmValue(iP, J, K) = mmA%l3mmValue(iP, J, K) + newField 
		     lonStart = lonStart - 2.0*PI
		     IF( lonStart < lons_rev(anlats(J, iP)) ) EXIT
		  ENDDO
		  IF( nc > 0) THEN
		     mmA%l3mmValue(iP, J, K) = mmA%l3mmValue(iP, J, K)/nc
		     mmA%l3mmPrecision(iP, J, K) = 0.0 
		  ELSE
		     mmA%l3mmValue(iP, J, K) = 0.0 
		     mmA%l3mmPrecision(iP, J, K) = 0.0 
		  END IF
		ENDDO

!*** Deallocate intermidiate arrays
                DEALLOCATE(Y2)
        	DEALLOCATE(lons_rev, fields_rev)

	     ELSE

                DO K = 1, cfProd%nLons
		  mmA%l3mmValue(iP, J, K) = 0.0 
		  mmA%l3mmPrecision(iP, J, K) = 0.0 
		ENDDO

	     END IF

	     ! Descending Mode

             IF( dnlats(J, iP) > 0 ) THEN

!*** Allocate space for all the arrays
        	ALLOCATE(lons_rev(dnlats(J, iP)), STAT=error)
        	ALLOCATE(fields_rev(dnlats(J, iP)), STAT=error)

                DO i = 1, dnlats(J, iP)
                    dlons(J, i, iP) = dlons(J, i, iP) + (i-1)*2.0*PI
                    lons_rev(i) = dlons(J, dnlats(J, iP)+1-i, iP)
                    fields_rev(i) = dfields(J, dnlats(J, iP)+1-i, iP)
                ENDDO 
                ALLOCATE(Y2(dnlats(J, iP)), STAT=error)
		CALL SPLINE(lons_rev, fields_rev, dnlats(J, iP), 0.0, 0.0, Y2)
                DO K = 1, cfProd%nLons
		  IF( cfProd%longGrid(K) > lons_rev(1) ) THEN
		     lonStart = cfProd%longGrid(K) - 2.0*PI
		  ELSE
		     lonStart = cfProd%longGrid(K)
		  END IF
		  nc = 0
		  mmD%l3mmValue(iP, J, K) = 0.0 
		  DO
		     nc = nc + 1
		     CALL SPLINT(lons_rev, fields_rev, Y2, dnlats(J, iP), lonStart, newField)
		     mmD%l3mmValue(iP, J, K) = mmD%l3mmValue(iP, J, K) + newField 
		     lonStart = lonStart - 2.0*PI
		     IF( lonStart < lons_rev(dnlats(J, iP)) ) EXIT
		  ENDDO
		  IF( nc > 0) THEN
		     mmD%l3mmValue(iP, J, K) = mmD%l3mmValue(iP, J, K)/nc
		     mmD%l3mmPrecision(iP, J, K) = 0.0 
		  ELSE
		     mmD%l3mmValue(iP, J, K) = 0.0 
		     mmD%l3mmPrecision(iP, J, K) = 0.0 
		  END IF
		ENDDO

!*** Deallocate intermidiate arrays
                DEALLOCATE(Y2)
        	DEALLOCATE(lons_rev, fields_rev)

	     ELSE

                DO K = 1, cfProd%nLons
		  mmD%l3mmValue(iP, J, K) = 0.0 
		  mmD%l3mmPrecision(iP, J, K) = 0.0 
		ENDDO

	     END IF

	     ! Combined Mode

             IF( anlats(J, iP) > 0 .or. dnlats(J, iP) > 0 ) THEN

                DO K = 1, cfProd%nLons
		  IF( anlats(J, iP) > 0 .and. dnlats(J, iP) == 0 ) THEN
		  	l3mm%l3mmValue(iP, J, K) = mmA%l3mmValue(iP, J, K)
		  	l3mm%l3mmPrecision(iP, J, K) = 0.0 
		  ELSE IF( anlats(J, iP) == 0 .and. dnlats(J, iP) > 0 ) THEN
		  	l3mm%l3mmValue(iP, J, K) = mmD%l3mmValue(iP, J, K)
		  	l3mm%l3mmPrecision(iP, J, K) = 0.0 
		  ELSE IF( anlats(J, iP) > 0 .and. dnlats(J, iP) > 0 ) THEN
		  	l3mm%l3mmValue(iP, J, K) = (mmA%l3mmValue(iP, J, K)+mmD%l3mmValue(iP, J, K))*0.5
		  	l3mm%l3mmPrecision(iP, J, K) = 0.0 
		  END IF
		ENDDO

	     ELSE

                DO K = 1, cfProd%nLons
		  l3mm%l3mmValue(iP, J, K) = 0.0 
		  l3mm%l3mmPrecision(iP, J, K) = 0.0 
		ENDDO

	     END IF

	  ENDDO

	ENDDO


!-----------------------------------
   END SUBROUTINE MonthlyMapFromL2
!-----------------------------------

!-------------------------------------------------------------------------
   INTEGER FUNCTION FindIndexForNormGrid(cfDef, alat)
!-------------------------------------------------------------------------

	Real (r8) alat
	TYPE( L3CFMDef_T ) :: cfDef

	integer i, j

	DO i = 1, cfDef%nNom 
	  IF(alat <= cfDef%l2nomLats(1)) THEN
		 FindIndexForNormGrid = 1
		 EXIT
	  ELSE IF(alat >= cfDef%l2nomLats(cfDef%nNom)) THEN
		 FindIndexForNormGrid = cfDef%nNom
		 EXIT
	  ELSE
		IF( i == cfDef%nNom ) THEN
			FindIndexForNormGrid = i
		 	EXIT
		ELSE IF( i == 1 ) THEN
			IF( alat <  (cfDef%l2nomLats(i)+cfDef%l2nomLats(i+1))*0.5 ) THEN 
				FindIndexForNormGrid = i
		 		EXIT
			END IF
		ELSE IF( alat >= (cfDef%l2nomLats(i)+cfDef%l2nomLats(i-1))*0.5 .and. 	&
			 alat <  (cfDef%l2nomLats(i)+cfDef%l2nomLats(i+1))*0.5 ) THEN
			FindIndexForNormGrid = i
		 	EXIT
		END IF
	  END IF	
	END DO

!-----------------------------------
   END FUNCTION FindIndexForNormGrid
!-----------------------------------

!-------------------------------------------------------------------------
   INTEGER FUNCTION FindLatIndexForL3Grid(cfProd, alat)
!-------------------------------------------------------------------------

	Real (r8) alat
        TYPE( L3CFMProd_T ) :: cfProd

	integer i, j

	DO i = 1, cfProd%nLats 
	  IF(alat <= cfProd%latGridMap(1)) THEN
		 FindLatIndexForL3Grid = 1
		 EXIT
	  ELSE IF(alat >= cfProd%latGridMap(cfProd%nLats)) THEN
		 FindLatIndexForL3Grid = cfProd%nLats 
		 EXIT
	  ELSE
		IF( i == cfProd%nLats ) THEN
			FindLatIndexForL3Grid = i
		 	EXIT
		ELSE IF( i == 1 ) THEN
			IF( alat <  (cfProd%latGridMap(i)+cfProd%latGridMap(i+1))*0.5 ) THEN 
				FindLatIndexForL3Grid = i
		 		EXIT
			END IF
		ELSE IF( alat >= (cfProd%latGridMap(i)+cfProd%latGridMap(i-1))*0.5 .and. 	&
			 alat <  (cfProd%latGridMap(i)+cfProd%latGridMap(i+1))*0.5 ) THEN
			FindLatIndexForL3Grid = i
		 	EXIT
		END IF
	  END IF	
	END DO

!-----------------------------------
   END FUNCTION FindLatIndexForL3Grid
!-----------------------------------

!-------------------------------------------------------------------------
   INTEGER FUNCTION FindLonIndexForL3Grid(cfProd, alon)
!-------------------------------------------------------------------------

	Real (r8) alon
        TYPE( L3CFMProd_T ) :: cfProd

	integer i, j

	DO i = 1, cfProd%nLons 
	  IF(alon <= cfProd%longGrid(1)) THEN
		 FindLonIndexForL3Grid = 1
		 EXIT
	  ELSE IF(alon >= cfProd%longGrid(cfProd%nLons)) THEN
		 FindLonIndexForL3Grid = cfProd%nLons 
		 EXIT
	  ELSE
		IF( i == cfProd%nLons ) THEN
			FindLonIndexForL3Grid = i
		 	EXIT
		ELSE IF( i == 1 ) THEN
			IF( alon <  (cfProd%longGrid(i)+cfProd%longGrid(i+1))*0.5 ) THEN 
				FindLonIndexForL3Grid = i
		 		EXIT
			END IF
		ELSE IF( alon >= (cfProd%longGrid(i)+cfProd%longGrid(i-1))*0.5 .and. 	&
			 alon <  (cfProd%longGrid(i)+cfProd%longGrid(i+1))*0.5 ) THEN
			FindLonIndexForL3Grid = i
		 	EXIT
		END IF
	  END IF	
	END DO

!-----------------------------------
   END FUNCTION FindLonIndexForL3Grid
!-----------------------------------

!-------------------------------------------------------------------------
   REAL FUNCTION FindRealLon(alon)
!-------------------------------------------------------------------------

	Real alon

	IF(mod(alon, 2.0*PI) < -PI) THEN
	   FindRealLon = 2.0*PI + mod(alon, 2.0*PI) 
	ELSE IF(mod(alon, 2.0*PI) > PI) THEN
	   FindRealLon = -2.0*PI + mod(alon, 2.0*PI) 
	ELSE
	   FindRealLon = mod(alon, 2.0*PI) 
	END IF
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


!-------------------------------------------------------------------------
   SUBROUTINE CalSolarRange(cfProd, l2gp, pStartIndex, pEndIndex,  	&
			    	anlats, dnlats, 		&
				asolarTime, dsolarTime, asolarZenith, dsolarZenith, 	&
			     	l3mz, mzA, mzD )
!-------------------------------------------------------------------------

! Brief description of program
! This is the main program to run the Core processing.

! Parameters

	! Variable definitions

        TYPE( PCFMData_T ) :: pcf
        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T ) :: cfDef
        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	TYPE( L3MZData_T ) :: l3mz, mzA, mzD

	Real, POINTER, DIMENSION(:, :, :) ::  asolarTime(:, :, :), dsolarTime(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  asolarZenith(:, :, :), dsolarZenith(:, :, :)

	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 

	Integer l2Days, pStartIndex, pEndIndex, iMin
	Integer iP, kP, J, I

	Real solarTimeMin, solarTimeMax, solarZenithMin, solarZenithMax

!*** Start calculation 

        iP = 0
	DO kP = pStartIndex, pEndIndex 
          iP = iP + 1

          DO J = 1, cfProd%nLats

	     ! Ascending Mode

             IF( anlats(J, iP) > 0 ) THEN

		solarTimeMin =  1.e20
		solarTimeMax = -1.e20
		solarZenithMin =  1.e20
		solarZenithMax = -1.e20

                DO I = 1, anlats(J, iP)
                    if(asolarTime(J, I, iP) >= solarTimeMax) then
			solarTimeMax = asolarTime(J, I, iP)
		    else if(asolarTime(J, I, iP) < solarTimeMin) then
			solarTimeMin = asolarTime(J, I, iP)
		    end if

                    if(asolarZenith(J, I, iP) >= solarZenithMax) then
			solarZenithMax = asolarZenith(J, I, iP)
		    else if(asolarZenith(J, I, iP) < solarZenithMin) then
			solarZenithMin = asolarZenith(J, I, iP)
		    end if
                ENDDO 
	     END IF

	     ! Descending Mode

             IF( dnlats(J, iP) > 0 ) THEN

		solarTimeMin =  1.e20
		solarTimeMax = -1.e20
		solarZenithMin =  1.e20
		solarZenithMax = -1.e20

                DO I = 1, dnlats(J, iP)
                    if(dsolarTime(J, I, iP) >= solarTimeMax) then
			solarTimeMax = dsolarTime(J, I, iP)
		    else if(dsolarTime(J, I, iP) < solarTimeMin) then
			solarTimeMin = asolarTime(J, I, iP)
		    end if

                    if(dsolarZenith(J, I, iP) >= solarZenithMax) then
			solarZenithMax = dsolarZenith(J, I, iP)
		    else if(dsolarZenith(J, I, iP) < solarZenithMin) then
			solarZenithMin = dsolarZenith(J, I, iP)
		    end if
                ENDDO 
	     END IF

          ENDDO 

        ENDDO 
!-----------------------------------
   END SUBROUTINE CalSolarRange
!-----------------------------------



!===================
END MODULE MonthlyProcessModule
!===================



