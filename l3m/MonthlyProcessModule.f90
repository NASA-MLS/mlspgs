
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
					l3mz, mzA, mzD, l3dz, dzA, dzD)
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
	Real, POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)

	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 

        CHARACTER (LEN=480) :: msr

	INTEGER ::  error, l2Days, nlev, nf, nwv, numDays, pEndIndex, pStartIndex

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
           dzA(j)%pressure = 0.0
           dzA(j)%latitude = 0.0
           dzA(j)%l3dzValue = 0.0
           dzA(j)%l3dzPrecision = 0.0
           dzD(j)%pressure = 0.0
           dzD(j)%latitude = 0.0
           dzD(j)%l3dzValue = 0.0
           dzD(j)%l3dzPrecision = 0.0

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

        mzA%pressure = 0.0
        mzA%latitude = 0.0
        mzA%l3mzValue = 0.0
        mzA%l3mzPrecision = 0.0

        mzD%pressure = 0.0
        mzD%latitude = 0.0
        mzD%l3mzValue = 0.0
        mzD%l3mzPrecision = 0.0

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

        mmA%pressure = 0.0
        mmA%latitude = 0.0
        mmA%longitude = 0.0
        mmA%l3mmValue = 0.0
        mmA%l3mmPrecision = 0.0

        mmD%pressure = 0.0
        mmD%latitude = 0.0
        mmD%longitude = 0.0
        mmD%l3mmValue = 0.0
        mmD%l3mmPrecision = 0.0

!*** Sort & Prepare the Data 

	Call SortDataMonthly(cfProd, l2Days, l2gp, 		&
			pStartIndex, pEndIndex,			&
			tau0, 					&
			anlats, dnlats, 			&
		        alats, dlats, 				&
			alons, dlons, 				&
			atimes, dtimes, 			&
			afields, dfields )

!*** Calculate Daily Zonal Means  

        Call DailyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days,        	&
                                    pStartIndex, pEndIndex,     &
                                    l3dz, dzA, dzD )

!*** Calculate Monthly Zonal Means  

        Call MonthlyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days,      	&
                                    pStartIndex, pEndIndex,     &
                                    l3mz, mzA, mzD )

!*** Calculate Monthly Map 

   	Call MonthlyMapFromL2(cfProd, l2gp, 		&
			    pStartIndex, pEndIndex,  	&
			    anlats, dnlats, 		&
		            alats, dlats, 		&
			    alons, dlons, 		&
			    atimes, dtimes, 		&
			    afields, dfields, 		&
			    l3mm, mmA, mmD )

!*** Allocate space for all the arrays

	DeAllocate(alats, dlats, alons, dlons, atimes, dtimes, afields, dfields)

!-----------------------------------
   END SUBROUTINE MonthlyCoreProcessing
!-----------------------------------



!-------------------------------------------------------------------------
   SUBROUTINE SortDataMonthly(cfProd, l2Days, l2gp, pStartIndex, pEndIndex, tau0, anlats, dnlats, 	&
		alats_interp, dlats_interp, alons_interp, dlons_interp, 	&
		atimes_interp, dtimes_interp, afields_interp, dfields_interp)
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex, dindex_prev

        TYPE( L3CFMProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	Real, POINTER, DIMENSION(:, :) ::  l2Times(:, :), l2Lons_new(:, :), l2Lons(:, :),  l2Lons_old(:, :), &
						l2Lats(:, :), l2Values(:, :)

	Real, POINTER, DIMENSION(:, :, :) ::  alons_interp(:, :, :), alats_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  dlons_interp(:, :, :), dlats_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  atimes_interp(:, :, :), dtimes_interp(:, :, :)
	Real, POINTER, DIMENSION(:, :, :) ::  afields_interp(:, :, :), dfields_interp(:, :, :)

	Real, POINTER, DIMENSION(:, :, :) ::  alons_interp_old(:, :, :), dlons_interp_old(:, :, :)

	Real dlons, alons, lons_found, lats_found, times_found, fields_found, 	&
		  slope, slope_time, slope_field, sTime, ddtad

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

        ALLOCATE(l2Lons_new(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Lons_old(nterms, pEndIndex-pStartIndex+1), STAT=error)

        ALLOCATE(alons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(alats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(atimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(afields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dtimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dfields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

        ALLOCATE(alons_interp_old(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlons_interp_old(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

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
		IF(ABS(slope) < 1.e-6) THEN
		   lons_found = l2Lons_new(lindex, iP)
		   lons_found_old = l2Lons_old(lindex, iP)
		ELSE
		   lons_found = l2Lons_new(lindex, iP) + (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		   lons_found_old = l2Lons_old(lindex, iP) + (cfProd%latGridMap(J)-l2Lats(lindex, iP))/slope
		END IF
		lats_found = cfProd%latGridMap(J)
		times_found = l2Times(lindex, iP) + (lons_found-l2Lons_new(lindex, iP))*slope_time
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
	           afields_interp(aindex, anlats(aindex, iP), iP) = fields_found
	           alons_interp_old(aindex, anlats(aindex, iP), iP) = lons_found_old
	         !END IF
	     END IF

	     IF( slope >= 0 .and. aindex /= dindex_prev) THEN
	         dindex_prev = aindex
		 IF(anlats(aindex, iP) > 0) THEN
       	           dnlats(aindex, iP) = dnlats(aindex, iP) + 1 
	           dlons_interp(aindex, dnlats(aindex, iP), iP) = lons_found
	           dlats_interp(aindex, dnlats(aindex, iP), iP) = lats_found
	           dtimes_interp(aindex, dnlats(aindex, iP), iP) = times_found
	           dfields_interp(aindex, dnlats(aindex, iP), iP) = fields_found
	           dlons_interp_old(aindex, dnlats(aindex, iP), iP) = lons_found_old
	         END IF
	     END IF


	    END IF

          ENDDO


        ENDDO

!*** Deallocate intermidiate arrays 

        DEALLOCATE(l2Lons, l2Lons_new, l2Lons_old, l2Lats, l2Values)

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

	Integer l2Days, pStartIndex, pEndIndex, iMin

	Real slope

!*** Allocate space for all the arrays

        ALLOCATE(iComArr(cfDef%nNom), STAT=error)
        ALLOCATE(iAscArr(cfDef%nNom), STAT=error)
        ALLOCATE(iDesArr(cfDef%nNom), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

	DO iT = 1, cfDef%nNom 
	  iComArr(iT) = 0
	  iAscArr(iT) = 0
	  iDesArr(iT) = 0
        ENDDO

!*** Re-arrange the data into longitude order for each pressure level 

        DO iD = 1, l2Days
	  iP = 0
	  DO kP = pStartIndex, pEndIndex 
	    iP = iP + 1
	    DO iT = 1, cfDef%nNom 
	      iComArr(iT) = 0
	      iAscArr(iT) = 0
	      iDesArr(iT) = 0

	      l3dz(iD)%l3dzValue(iP, iT) = 0.0 
	      dzA(iD)%l3dzValue(iP, iT) = 0.0 
	      dzD(iD)%l3dzValue(iP, iT) = 0.0 
            ENDDO
	    DO iT = 1, l2gp(iD)%nTimes-1
	      !iCom = real(l2gp(iD)%latitude(iT)-cfDef%l2nomLats(1))/real(cfDef%l2nomLats(2)-cfDef%l2nomLats(1))+1.5
   	      iCom = FindIndexForNormGrid(cfDef, l2gp(iD)%latitude(iT))
	      l3dz(iD)%l3dzValue(iP, iCom) = l3dz(iD)%l3dzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	      iComArr(iCom) = iComArr(iCom) + 1 
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
	      ELSE
	        dzD(iD)%l3dzValue(iP, iCom) = dzD(iD)%l3dzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iDesArr(iCom) = iDesArr(iCom) + 1 
	      END IF	
	    ENDDO

	    DO iT = 1, cfDef%nNom 
	      IF(iComArr(iT) == 0) THEN
	      	l3dz(iD)%l3dzValue(iP, iT) = 0.0 
	      	dzA(iD)%l3dzValue(iP, iT) = 0.0 
	      	dzD(iD)%l3dzValue(iP, iT) = 0.0 
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
	      END IF
	    ENDDO

	  ENDDO
	ENDDO

!*** Deallocate 

        DEALLOCATE(iComArr, iAscArr, iDesArr)

!-----------------------------------
   END SUBROUTINE DailyZonalMeanFromL2
!-----------------------------------

!-------------------------------------------------------------------------
   SUBROUTINE MonthlyZonalMeanFromL2(cfProd, cfDef, l2gp, l2Days, pStartIndex, pEndIndex,  	&
			     	l3mz, mzA, mzD )
!-------------------------------------------------------------------------

        integer error, i, j, iT, iD, iP, kP, nterms, nstart, nr
	INTEGER ::  iCom, iAasc, iDes, nloop, aindex

        TYPE( L3CFMProd_T ) :: cfProd
	TYPE( L3CFMDef_T ) :: cfDef
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
	TYPE( L3MZData_T ) :: l3mz, mzA, mzD

	Integer, POINTER, DIMENSION(:) ::  iComArr(:), iAscArr(:), iDesArr(:)

	Integer l2Days, pStartIndex, pEndIndex, iMin

	Real slope

!*** Allocate space for all the arrays

        ALLOCATE(iComArr(cfDef%nNom), STAT=error)
        ALLOCATE(iAscArr(cfDef%nNom), STAT=error)
        ALLOCATE(iDesArr(cfDef%nNom), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

	DO iT = 1, cfDef%nNom 
	  iComArr(iT) = 0
	  iAscArr(iT) = 0
	  iDesArr(iT) = 0
        ENDDO

!*** Re-arrange the data into longitude order for each pressure level 

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
          DO iD = 1, l2Days
	    DO iT = 1, l2gp(iD)%nTimes-1
	      !iCom = real(l2gp(iD)%latitude(iT)-l3mz%latitude(1))/real(l3mz%latitude(2)-l3mz%latitude(1))+1.5
   	      iCom = FindIndexForNormGrid(cfDef, l2gp(iD)%latitude(iT))
	      l3mz%l3mzValue(iP, iCom) = l3mz%l3mzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	      iComArr(iCom) = iComArr(iCom) + 1 
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
	      ELSE
	        mzD%l3mzValue(iP, iCom) = mzD%l3mzValue(iP, iCom) + l2gp(iD)%l2gpValue(1, kP, iT)
	        iDesArr(iCom) = iDesArr(iCom) + 1 
	      END IF	
	    ENDDO
	  ENDDO

	  DO iT = 1, cfDef%nNom 
	    IF(iComArr(iT) == 0) THEN
	      	l3mz%l3mzValue(iP, iT) = 0.0 
	      	mzA%l3mzValue(iP, iT) = 0.0 
	      	mzD%l3mzValue(iP, iT) = 0.0 
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
	    END IF
	  ENDDO

	ENDDO

!*** Deallocate 

        DEALLOCATE(iComArr, iAscArr, iDesArr)

!-----------------------------------
   END SUBROUTINE MonthlyZonalMeanFromL2
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

        
!*** Start calculation 

        iP = 0
        DO kP = pStartIndex, pStartIndex
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
		ELSE IF( alat > cfDef%l2nomLats(i) .and. alat < cfDef%l2nomLats(i+1) ) THEN
			FindIndexForNormGrid = i
		 	EXIT
		END IF
	  END IF	
	END DO

!-----------------------------------
   END FUNCTION FindIndexForNormGrid
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



!===================
END MODULE MonthlyProcessModule
!===================



