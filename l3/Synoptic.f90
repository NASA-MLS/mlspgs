
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
   USE L3DZData
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
   SUBROUTINE DailyCoreProcessing(cfProd, pcf, l2Days, l2gp, avgPeriod, l3sp, l3dm, dmA, dmD, &
					l3dz, dzA, dzD, l3r, residA, residD, flags)
!-------------------------------------------------------------------------

! Brief description of program
! This is the main program to run the Core processing.

! Parameters

	integer, Parameter :: nlons = 180

	! Variable definitions

        TYPE( PCFData_T ) :: pcf
        TYPE( L3CFProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)
        TYPE( L3DMData_T ), POINTER :: l3dm(:), dmA(:), dmD(:)
	TYPE( L3DZData_T ), POINTER :: l3dz(:), dzA(:), dzD(:)
        TYPE( L3SPData_T ), POINTER :: l3sp(:)
        TYPE( L2GPData_T ), POINTER :: l3r(:), residA(:), residD(:)

   	TYPE( OutputFlags_T ) :: flags

	Real (r8), POINTER, DIMENSION(:, :, :) ::  alons(:, :, :), alats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  dlons(:, :, :), dlats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  atimes(:, :, :), dtimes(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)

	Real (r8), POINTER, DIMENSION(:) :: l2gpValue_1Day
        REAL(r8), POINTER :: avgPeriod(:)

	INTEGER, DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: anlats, dnlats 
	Real (r8), DIMENSION(cfProd%nLats, l2gp(1)%nLevels) :: delTad 
        REAL(r8), DIMENSION(:), POINTER :: l3Result       ! returned reconstructed result 


        CHARACTER (LEN=480) :: msr

	INTEGER ::  error, l2Days, nlev, nf, nwv, numDays, numSwaths, rDays, pEndIndex, pStartIndex

        integer i, j, k, ctype, iP, kP, iD, iL
        real(r8), Dimension(nlons) :: xlon, result
	real zonAvg, tau0, l3ret

!*** Initilize variables
 
	nlev = 24
        nwv = 10
        nf = cfProd%nFreqs

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

        IF (cfProd%mode == 'all') THEN
           numSwaths = 3
        ELSE
           numSwaths = 1
        ENDIF

        ALLOCATE( l3sp(numSwaths), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3sp array.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

        IF (cfProd%mode == 'all') THEN
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
           l3sp(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           l3sp(j)%latitude = cfProd%latGridMap(:l3sp(j)%nLats)

        ENDDO

        numDays = cfProd%nDays
        ALLOCATE( l3dm(numDays), dmA(numDays), dmD(numDays), STAT=error )
        ALLOCATE( l3dz(numDays), dzA(numDays), dzD(numDays), STAT=error )
        IF ( error /= 0 ) THEN
           msr = MLSMSG_Allocate // ' l3dm, l3r arrays.'
           CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF

!!      Initialize Daily Map 

        l3dm%name = cfProd%l3prodNameD
        dmA%name = TRIM(cfProd%l3prodNameD) // 'Ascending'
        dmD%name = TRIM(cfProd%l3prodNameD) // 'Descending'

        DO j = 1, numDays

           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, l3dm(j) )
           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, dmA(j) )
           CALL AllocateL3DM( nlev, cfProd%nLats, cfProd%nLons, dmD(j) )

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

        ENDDO

!!      Initialize Daily Zonal Mean 

        l3dz%name = cfProd%l3prodNameD
        dzA%name = TRIM(cfProd%l3prodNameD) // 'Ascending'
        dzD%name = TRIM(cfProd%l3prodNameD) // 'Descending'

        DO j = 1, numDays

           CALL AllocateL3DZ( nlev, cfProd%nLats, l3dz(j) )
           CALL AllocateL3DZ( nlev, cfProd%nLats, dzA(j) )
           CALL AllocateL3DZ( nlev, cfProd%nLats, dzD(j) )

           l3dz(j)%time = cfProd%timeD(j)
           dzA(j)%time = cfProd%timeD(j)
           dzD(j)%time = cfProd%timeD(j)

           l3dz(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           dzA(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)
           dzD(j)%pressure = l2gp(1)%pressures(pStartIndex:pEndIndex)

           l3dz(j)%latitude = cfProd%latGridMap(:l3dm(j)%nLats)
           dzA(j)%latitude = cfProd%latGridMap(:dmA(j)%nLats)
           dzD(j)%latitude = cfProd%latGridMap(:dmD(j)%nLats)

        ENDDO

!!      Initialize Daily Map Residues

        CALL ReadL2GPProd(cfProd, pcf%l3StartDay, pcf%l3EndDay, rDays, l3r)
        CALL ReadL2GPProd(cfProd, pcf%l3StartDay, pcf%l3EndDay, rDays, residA)
        CALL ReadL2GPProd(cfProd, pcf%l3StartDay, pcf%l3EndDay, rDays, residD)

        l3r%name    = TRIM(cfProd%l3prodNameD) // 'Residuals'
        residA%name = TRIM(cfProd%l3prodNameD) // 'AscendingResiduals'
        residD%name = TRIM(cfProd%l3prodNameD) // 'DescendingResiduals'

        DO j = 1, rDays
           l3r(j)%l2gpValue    = 0.0
           residA(j)%l2gpValue = 0.0
           residD(j)%l2gpValue = 0.0
        ENDDO
 
!!      Initialize Flags

        flags%writel3sp    = .FALSE.
        flags%writel3dmCom = .FALSE.
        flags%writel3dmAsc = .FALSE.
        flags%writel3dzAsc = .FALSE.
        flags%writel3dmDes = .FALSE.
	flags%writel3dzDes = .FALSE.
        flags%writel3rCom  = .FALSE.
        flags%writel3rAsc  = .FALSE.
        flags%writel3rDes  = .FALSE.
 

!*** Calculate average orbital period (day)

	tau0 = 0.0
        DO I = 1, size(avgPeriod)
	  tau0 = tau0 + avgPeriod(i)
	ENDDO 
	tau0 = tau0/86400.0/float(size(avgPeriod))
 
!*** Sort & Prepare the Data 

	Call SortData(cfProd, pcf, l2Days, l2gp, 		&
			pStartIndex, pEndIndex,			&
			tau0, 					&
			anlats, dnlats, 			&
		        alats, dlats, 				&
			alons, dlons, 				&
			atimes, dtimes, 			&
			afields, dfields, delTad )

	!open(5, file='l2gp_intp.dat', status='replace')
	!DO J = 1, anlats(10, 1)
	!  write(5,*) J, atimes(10, J, 1)*60.*60.*24.+l2gp(1)%time(1), alons(10, J, 1), afields(10, J, 1) 
	  !write(5,*) J, atimes(10, J, 1), alons(10, J, 1), afields(10, J, 1) 
        !ENDDO
        !close(5)

!*** Main Loop 

!	cfProd%mode = 'com'

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO J = 1, cfProd%nLats
	  !DO J = 10, 10 
       	     IF( anlats(J, iP) > 0 ) THEN
	        CALL Init(cfProd%mode, 				&
			  nt_a_i   = anlats(J, iP), 		&
			  nt_d_i   = dnlats(J, iP), 		&
			  tau0_i   = tau0, 			&
			  delTad_i = delTad(J, iP), 		&
			  c0_i     = 2.0*PI, 			&
			  lonD0_i  = alons(J, 1, iP), 		&
			  tD0_i    = atimes(J, 1, iP), 		&
			  lonA0_i  = dlons(J, 1, iP), 		&
			  tA0_i    = dtimes(J, 1, iP), 		&
			  lat_i    = alats(J, 1, iP) )
          	CALL DataGenerate(cfProd%mode, afields(J, :, iP), dfields(J, :, iP) )
		IF (cfProd%mode == 'com') THEN
                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      l3dz(iD)%l3dzValue(iP, J) = zonAvg/float(l3dm(iD)%nLons) 
		   ENDDO
                   DO iL = 1, anlats(J, iP)
		      !CALL Diagnostics(cfProd%mode, atimes(J, iL, iP), alons(J, iL, iP), l3ret) 
		      !print *, atimes(J, iL, iP), alons(J, iL, iP), afields(J, iL, iP)-l3ret
		   ENDDO
                   DO iL = 1, anlats(J, iP)
		      !CALL Diagnostics(cfProd%mode, dtimes(J, iL, iP), dlons(J, iL, iP), l3ret) 
		      !print *, dfields(J, iL, iP)-l3ret
		   ENDDO
		   DeAllocate(l3Result)
		   flags%writel3dmCom = .TRUE.
		   flags%writel3dzCom = .TRUE.
		   flags%writel3sp = .TRUE.
		ELSE IF (cfProd%mode == 'asc') THEN
                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      dzA(iD)%l3dzValue(iP, J) = zonAvg/float(dmA(iD)%nLons) 
		   ENDDO
		   DeAllocate(l3Result)
		   flags%writel3dmAsc = .TRUE.
		   flags%writel3dzAsc = .TRUE.
		   flags%writel3sp = .TRUE.
		ELSE IF (cfProd%mode == 'des') THEN
                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct(cfProd%mode, real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      dzD(iD)%l3dzValue(iP, J) = zonAvg/float(dmD(iD)%nLons) 
		   ENDDO
		   DeAllocate(l3Result)
		   flags%writel3dmDes = .TRUE.
		   flags%writel3dzDes = .TRUE.
		   flags%writel3sp = .TRUE.
		ELSE IF (cfProd%mode == 'all') THEN
                   ALLOCATE(l3Result(l3dm(1)%nLons), STAT=error)
                   CALL CordTransform('com')
	   	   CALL FFSM(l3sp(1), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('com', real(l3dm(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    l3dm(iD)%nLons, l3dm(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, l3dm(iD)%nLons
			l3dm(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      l3dz(iD)%l3dzValue(iP, J) = zonAvg/float(l3dm(iD)%nLons) 
		   ENDDO
		   DeAllocate(l3Result)

		   flags%writel3dmCom = .TRUE.
		   flags%writel3dzCom = .TRUE.

                   ALLOCATE(l3Result(dmA(1)%nLons), STAT=error)
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMA(l3sp(2), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('asc', real(dmA(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmA(iD)%nLons, dmA(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, dmA(iD)%nLons
			dmA(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      dzA(iD)%l3dzValue(iP, J) = zonAvg/float(dmA(iD)%nLons) 
		   ENDDO
		   DeAllocate(l3Result)

		   flags%writel3dmAsc = .TRUE.
		   flags%writel3dzAsc = .TRUE.

                   ALLOCATE(l3Result(dmD(1)%nLons), STAT=error)
                   CALL CordTransform(cfProd%mode)
	   	   CALL FFSMD(l3sp(3), iP, J)
                   DO iD = 1, cfProd%nDays
        	      CALL Reconstruct('des', real(dmD(iD)%time-l2gp(1)%time(1))/86400.0, 	&
				    dmD(iD)%nLons, dmD(iD)%longitude, l3Result)
		      zonAvg = 0.0
                      DO I = 1, dmD(iD)%nLons
			dmD(iD)%l3dmValue(iP, J, I) = l3Result(I) 
		        zonAvg = zonAvg + l3Result(I) 
		      ENDDO
		      dzD(iD)%l3dzValue(iP, J) = zonAvg/float(dmD(iD)%nLons) 
		   ENDDO
		   DeAllocate(l3Result)

		   flags%writel3dmDes = .TRUE.
		   flags%writel3dzDes = .TRUE.
		   flags%writel3sp = .TRUE.
		END IF

                CALL ClearMemory()

	     END IF
	  ENDDO

	ENDDO

!-----------------------------------
   END SUBROUTINE DailyCoreProcessing
!-----------------------------------


!-------------------------------------------------------------------------
   SUBROUTINE ParameterCalc(cfProd, pcf, l2Days, l2gp)
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, ctype, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex, dindex_prev

        TYPE( PCFData_T ) :: pcf
        TYPE( L3CFProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)


!-----------------------------------
   END SUBROUTINE ParameterCalc
!-----------------------------------


!-------------------------------------------------------------------------
   SUBROUTINE SortData(cfProd, pcf, l2Days, l2gp, pStartIndex, pEndIndex, tau0, anlats, dnlats, 	&
		alats_interp, dlats_interp, alons_interp, dlons_interp, 	&
		atimes_interp, dtimes_interp, afields_interp, dfields_interp,	&
		delTad)
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, kP, ctype, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex, dindex_prev

        TYPE( PCFData_T ) :: pcf
        TYPE( L3CFProd_T ) :: cfProd
        TYPE( L2GPData_T ), POINTER :: l2gp(:)

	Real (r8), POINTER, DIMENSION(:, :) ::  l2Times(:, :), l2Times_Rev(:, :), &
						l2Lons_new(:, :), l2Lons(:, :), l2Lons_Rev(:, :), &
						l2Lats(:, :), l2Lats_Rev(:, :), &
						l2Values(:, :), l2Values_Rev(:, :)

	Real (r8), POINTER, DIMENSION(:, :, :) ::  alons_interp(:, :, :), alats_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  dlons_interp(:, :, :), dlats_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  atimes_interp(:, :, :), dtimes_interp(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  afields_interp(:, :, :), dfields_interp(:, :, :)

	Real (r8) dlons, alons, lons_found, lats_found, times_found, fields_found, 	&
		  slope, slope_time, slope_field, sTime, ddtad

	Real tau0
 
	Integer nPd, pStartIndex, pEndIndex

	INTEGER, DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: anlats, dnlats
	Real (r8), DIMENSION(cfProd%nLats, pEndIndex-pStartIndex+1) :: delTad 

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
        ALLOCATE(l2Lons_Rev(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Times_Rev(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Lats_Rev(nterms, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(l2Values_Rev(nterms, pEndIndex-pStartIndex+1), STAT=error)

        ALLOCATE(alons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(alats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(atimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(afields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlons_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dlats_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dtimes_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)
        ALLOCATE(dfields_interp(cfProd%nLats, nPd*l2Days, pEndIndex-pStartIndex+1), STAT=error)

	IF(error /= 0 ) THEN
          print *, "Allocation Error"
	  STOP
	END IF

!*** Re-arrange the data into longitude order for each pressure level 

!	open(5, file='l2gp_info.dat', status='replace')

        nstart = 1
	DO iD = 1, l2Days
	  DO iT = 1, l2gp(iD)%nTimes
	    iP = 0
	    DO kP = pStartIndex, pEndIndex 
	        iP = iP + 1
                l2Times(nstart, iP) = l2gp(iD)%time(iT)
                l2Lons(nstart, iP) = l2gp(iD)%longitude(iT) * PI/180.0
                l2Lats(nstart, iP) = l2gp(iD)%latitude(iT)
                l2Values(nstart, iP) = l2gp(iD)%l2gpValue(1, kP, iT) 
	    ENDDO
!	    write(5, *) nstart, l2Times(nstart, 1), l2Lons(nstart, 1), l2Lats(nstart, 1), l2Values(nstart, 1)
            nstart = nstart + 1
	  ENDDO
	ENDDO

!	close(5)
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

!*** Reverse longitude order

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO iT = 1, nterms 
     		l2Lons_Rev(iT, iP) = l2Lons_new(nterms+1-iT, iP)
	  ENDDO
	ENDDO

!*** Reverse time order

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
          sTime = l2Times(1, iP)
	  DO iT = 1, nterms 
                l2Times(iT, iP) = (l2Times(iT, iP) - sTime)/86400.0 
	  ENDDO
	ENDDO

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO iT = 1, nterms 
                l2Times_Rev(iT, iP) = l2Times(nterms+1-iT, iP)
	  ENDDO
	ENDDO

!*** Reverse latitude order

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO iT = 1, nterms 
     		l2Lats_Rev(iT, iP) = l2Lats(nterms+1-iT, iP)
	  ENDDO
	ENDDO

!*** Reverse field value order

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1
	  DO iT = 1, nterms 
     		l2Values_Rev(iT, iP) = l2Values(nterms+1-iT, iP)
	  ENDDO
	ENDDO

!*** Deallocate intermidiate arrays 

        DEALLOCATE(l2Lons, l2Lons_new, l2Lats, l2Values)

!*** Find longitudes & latitudes for L3 latitudes in L2 data points 

        dlons = 0.3*PI/180.0

	iP = 0
	DO kP = pStartIndex, pEndIndex 
	  iP = iP + 1

 	  alons = l2Lons_Rev(1, iP) 

	  lindex = 1
	  lindex_prev = 0
	  aindex_prev = 0
	  dindex_prev = 0

          nloop = (l2Lons_Rev(nterms, iP)-l2Lons_Rev(1, iP))/dlons

	  DO J = 1, cfProd%nLats
       	     anlats(J, iP) = 0 
       	     dnlats(J, iP) = 0 
	  ENDDO

	  DO I = 1, nloop
	    alons = alons +  dlons
	    DO iT = lindex, nterms-1 
              IF( l2Lons_Rev(iT, iP) < alons .AND. l2Lons_Rev(iT+1, iP) >= alons ) THEN
	 	 lindex = iT
	 	 EXIT
	      END IF
	    ENDDO

            IF ( lindex /= lindex_prev ) THEN

	     lindex_prev = lindex

	     DO J = 1, cfProd%nLats
              IF( ( l2Lats_Rev(lindex, iP) <= cfProd%latGridMap(J) .AND. 	&
		    l2Lats_Rev(lindex+1, iP) > cfProd%latGridMap(J) )  .OR.	&
                  ( l2Lats_Rev(lindex, iP) >= cfProd%latGridMap(J) .AND. 	&
		    l2Lats_Rev(lindex+1, iP) < cfProd%latGridMap(J) )  ) THEN
		slope = (l2Lats_Rev(lindex+1, iP)-l2Lats_Rev(lindex, iP))/	&
			(l2Lons_Rev(lindex+1, iP)-l2Lons_Rev(lindex, iP))
		slope_time = (l2Times_Rev(lindex+1, iP)-l2Times_Rev(lindex, iP))/	&
			(l2Lons_Rev(lindex+1, iP)-l2Lons_Rev(lindex, iP))
		slope_field = (l2Values_Rev(lindex+1, iP)-l2Values_Rev(lindex, iP))/	&
			(l2Lons_Rev(lindex+1, iP)-l2Lons_Rev(lindex, iP))
		IF(ABS(slope) < 1.e-6) THEN
		   lons_found = l2Lons_Rev(lindex, iP)
		ELSE
		   lons_found = l2Lons_Rev(lindex, iP) + (cfProd%latGridMap(J)-l2Lats_Rev(lindex, iP))/slope
		END IF
		lats_found = cfProd%latGridMap(J)
		times_found = l2Times_Rev(lindex, iP) + (lons_found-l2Lons_Rev(lindex, iP))*slope_time
		fields_found = l2Values_Rev(lindex, iP) + (lons_found-l2Lons_Rev(lindex, iP))*slope_field
                aindex = J
                if(aindex == 1) then
	           !print *, 'aindex=1 case:', lats_found, lons_found, times_found
	  	end if
		EXIT
	      END IF
	     ENDDO
              
             !if(lats_found <= -74.0 .and. lons_found < -2830.0) then
	!	print *, aindex, lons_found, lats_found
	!     end if

	     IF( slope <= 0 .and. aindex /= aindex_prev) THEN
	         aindex_prev = aindex
       	         anlats(aindex, iP) = anlats(aindex, iP) + 1 
	         alons_interp(aindex, anlats(aindex, iP), iP) = lons_found
	         alats_interp(aindex, anlats(aindex, iP), iP) = lats_found
	         atimes_interp(aindex, anlats(aindex, iP), iP) = times_found
	         afields_interp(aindex, anlats(aindex, iP), iP) = fields_found
	     END IF

	     IF( slope > 0 .and. aindex /= dindex_prev) THEN
	         dindex_prev = aindex
       	         dnlats(aindex, iP) = dnlats(aindex, iP) + 1 
	         dlons_interp(aindex, dnlats(aindex, iP), iP) = lons_found
	         dlats_interp(aindex, dnlats(aindex, iP), iP) = lats_found
	         dtimes_interp(aindex, dnlats(aindex, iP), iP) = times_found
	         dfields_interp(aindex, dnlats(aindex, iP), iP) = fields_found
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
	    DO I = 1, anlats(J, iP) 
                delTad(J,iP) = delTad(J,iP) + atimes_interp(J, I, iP) - dtimes_interp(J, I, iP)
	    ENDDO
            IF(anlats(J, iP) > 0) THEN
              delTad(J,iP) = delTad(J,iP)/float(anlats(J, iP))
	    ELSE
              delTad(J,iP) = 0.0 
	    END IF
            IF(delTad(J,iP) < 0) THEN
              delTad(J,iP) = delTad(J,iP) + tau0 
	    END IF
	  ENDDO
	ENDDO

!*** Deallocate intermidiate arrays 


!-----------------------------------
   END SUBROUTINE SortData
!-----------------------------------


!-------------------------------------------------------------------------
   SUBROUTINE TransformCord(cfProd, pcf, l2Days, alats, dlats, alons, dlons, atimes, dtimes, afields, dfields)
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, ctype, nterms, nstart, nr, lindex, lindex_prev
	INTEGER ::  l2Days, nloop, aindex, aindex_prev, dindex, dindex_prev

        TYPE( PCFData_T ) :: pcf
        TYPE( L3CFProd_T ) :: cfProd

	Real (r8), POINTER, DIMENSION(:, :) ::  l2Times(:, :), l2Times_Rev(:, :), &
						l2Lons_new(:, :), l2Lons(:, :), l2Lons_Rev(:, :), &
						l2Lats(:, :), l2Lats_Rev(:, :), &
						l2Values(:, :), l2Values_Rev(:, :)

	Real (r8), POINTER, DIMENSION(:, :, :) ::  alons(:, :, :), alats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  dlons(:, :, :), dlats(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  atimes(:, :, :), dtimes(:, :, :)
	Real (r8), POINTER, DIMENSION(:, :, :) ::  afields(:, :, :), dfields(:, :, :)


 
! Read the l2gp data for ClO, the first product listed in the cf (and the only
! one for which simulated input files currently exist)



!-----------------------------------
   END SUBROUTINE TransformCord
!-----------------------------------



!-------------------------------------------------------------------------
   INTEGER FUNCTION FindL2LongitudeIndex(lons, alons)
!-------------------------------------------------------------------------

        integer error, i, j, k, iT, iD, iP, ctype, nterms, nstart, nr

        Real (r8) alons

	Real (r8), POINTER, DIMENSION(:, :) ::  lons(:, :)

 	FindL2LongitudeIndex = 1	
!-----------------------------------
   END FUNCTION FindL2LongitudeIndex
!-----------------------------------



!===================
END MODULE Synoptic
!===================

! $Log$
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


