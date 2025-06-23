! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
MODULE L2Interface
!==============================================================================

  USE intrinsic, ONLY: l_swath
  USE L2GPData, ONLY: L2GPData_T, ReadL2GPData, WriteL2GPData, SetupNewL2GPRecord, &
       & DestroyL2GPContents
  USE MLSCommon, ONLY: r8, MLSFile_T
  USE MLSFiles, ONLY: InitializeMLSFile, mls_openFile, mls_closeFile, & 
       & mls_hdf_version, mls_inqswath, MLS_CloseFile, MLS_OpenFile
  USE MLSL3Common, ONLY: DATE_LEN, CCSDS_LEN, FILENAMELEN, maxWindow, & 
       & GRIDNAMELEN, maxMisDays
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
       & MLSMSG_WARNING, MLSMSG_INFO, MLSMSG_FILEOPEN
  USE MLSPCF3, ONLY: mlspcf_l2gp_start, mlspcf_l2gp_end, &
       & mlspcf_l3dm_start, mlspcf_l3dm_end
  USE output_m, only: output
  USE PCFModule, ONLY: FindFileDay, ExpandFileTemplate
  USE Dump_0, only: DUMP

  IMPLICIT NONE
  private
  PUBLIC :: GetL2GPfromPCF, ReadL2GP, ReadL2GPProd, ResidualOutput, &
	& ReadL2DGData, ReadL2GPAttribute, SetupL2GPProd
  
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  
  ! Contents:
  
  ! Subroutines -- GetL2GPfromPCF
  !                ReadL2GP
  !                ReadL2GPProd
  !                ResidualOutput
  !                ReadL2DGData
  !                ReadL2GPAttribute
  
  ! Remarks:  This module contains subroutines involving the interface between
  !           Level 3 and L2.
  
  ! Parameters
  
CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELD1 = 'L2gpValue'
CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELD2 = 'L2gpPrecision'
CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELD3 = 'Status'
CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELD4 = 'Quality'
CHARACTER (LEN=*), PARAMETER, PUBLIC :: DATA_FIELD5 = 'L3Residual'
  
CONTAINS
  
  !----------------------------------------------------------------------------
  SUBROUTINE GetL2GPfromPCF(mlspcf_start, mlspcf_end, template, startDOY, &
       & endDOY, numFiles, pcfNames, mis_numFiles, mis_Days)
  !----------------------------------------------------------------------------
  USE SDPToolkit, ONLY: PGS_S_SUCCESS, Pgs_pc_getReference
    
    ! Brief description of subroutine
    ! This routine searches the PCF for l2gp files for a desired species.  It
    ! returns an array of names and the number of the files found.
    
    ! Arguments
    
    CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY
    CHARACTER (LEN=*), INTENT(IN) :: template
    INTEGER, INTENT(IN) :: mlspcf_end, mlspcf_start
    INTEGER, INTENT(OUT) :: numFiles, mis_numFiles
    INTEGER, INTENT(OUT) :: mis_Days(:)
    CHARACTER (LEN=FileNameLen) :: pcfNames(:)
   
    ! Parameters
    
    ! Functions
    
    ! Variables
    
    CHARACTER (LEN=FileNameLen) :: physicalFilename, type
    CHARACTER (LEN=8) :: date
    CHARACTER (LEN=7) :: dateChar
    CHARACTER (LEN=3) :: startChar, endChar, dayChar
    
    INTEGER :: i, j, indx, returnStatus, version, err
    INTEGER :: startInt, endInt, dateInt, dayInt
    INTEGER :: doy_current, doy_prev, day_diff
    LOGICAL, parameter :: DEBUG = .false. 
  
    ! Expand the level/species template
    
    CALL ExpandFileTemplate(template, type, 'L2GP')

    if (DEBUG) then
       call output('type= ', advance='no')   
       call output(trim(type), advance='yes')   
    end if

    numFiles = 0
    mis_numFiles = 0
    mis_Days = 0
    doy_current = 0 
    doy_prev = 0 
    day_diff = 0 

    ! Loop through all the PCF numbers for L2GP files
   
    do i = mlspcf_start, mlspcf_end
       version = 1
       returnStatus = Pgs_pc_getReference(i, version, physicalFilename)

       if (i == mlspcf_start) then
           indx = INDEX(physicalFilename, '.', .TRUE.)
           read(physicalFilename(indx-3:indx-1),'(i3)') doy_current
	   doy_prev = doy_current
       endif 
       if (DEBUG) then
         call output( 'physicalFilename = ', advance='no')
         call output( trim(physicalFilename), advance='yes')
         call output( 'returnStatus = ', advance='no')
         call output( returnStatus, advance='yes')
       end if

       ! If no file name was returned, go on to the next PCF number
    
       if(returnStatus /= PGS_S_SUCCESS) THEN
          if ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN
                ! Extract the date from the file name
                indx = INDEX(physicalFilename, '.', .TRUE.)
                date = physicalFilename(indx-8:indx-5)//'-'//physicalFilename(indx-3:indx-1)
                ! Check that the date is within the desired boundaries for reading
             
                if ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN
	           startChar = startDOY(6:8)
                   read (startChar, '(i3)') startInt
	           endChar = endDOY(6:8)
                   read (endChar, '(i3)') endInt
		   
		   if (numFiles <= (endInt - startInt)) then 
               
                   ! Save the name in an array of files for the species
                       mis_numFiles = mis_numFiles + 1
                
		   ! convert from char to int. There was problem using charstring 
		   ! with tk 5.2.10 (hdfeos5.1.6.2) 3/2004                         
	              dateChar = date(1:4) // date(6:8) 
                      read (dateChar, '(i7)') dateInt
		      dateInt = dateInt + 1
                      mis_Days(mis_numFiles) = dateInt
		   endif
                endif
          endif
          CYCLE
       endif 
       
       ! Check that the returned file name is for the proper species
       if ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN
          ! Extract the date from the file name
          indx = INDEX(physicalFilename, '.', .TRUE.)
          date = physicalFilename(indx-8:indx-5)//'-'//physicalFilename(indx-3:indx-1)
          
          ! Check that the date is within the desired boundaries for reading
	  if (DEBUG) then
             call output('date=', advance='no')
	     call output( date, advance='yes')
	  end if

          if ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN
             ! Save the name in an array of files for the species
             numFiles = numFiles + 1
             pcfNames(numFiles) = physicalFilename
	     if (DEBUG) then
   	     	call output('numFiles = ', advance ='no')          
   	     	call output(numFiles, advance ='yes')          
   	     	call output('pcfNames = ', advance ='no')          
   	     	call output(pcfNames(numFiles), advance ='yes')          
	     endif
          endif
       
	  ! Check missing days that are not in the PCF

          read(physicalFilename(indx-3:indx-1),'(i3)') doy_current
	  if (DEBUG) then
             print *,' doy_current: ', doy_current
             print *,' filename here: ', trim(physicalFilename)
             print *,' type: ', trim(type)
             print *,' i: ',i 
	  endif
	  day_diff = doy_current - doy_prev
	  if (day_diff > 1) then
	     do j = 1, day_diff - 1 
		mis_numFiles = mis_numFiles + 1
	        mis_Days(mis_numFiles) = doy_prev + j
	     enddo
	  endif
          doy_prev = doy_current
       endif
    enddo
   
  !-------------------------------
  END SUBROUTINE GetL2GPfromPCF
  !-------------------------------
  
  !--------------------------------------------------------
  SUBROUTINE ReadL2GPAttribute(l3StartDay, l3EndDay, tempfile, nDays, L3MFlag)
  !--------------------------------------------------------
    
    ! Brief description of subroutine
    ! This routine reads an L2GP file, returning a filled data structure.

    USE HDF, only: DFACC_READ 
    !USE MLSHDFEOS, only: he5_EHrdglatt
    USE SDPToolkit, only: max_orbits
    USE Allocate_Deallocate, only: ALLOCATE_TEST
    USE PCFHdr, only: GlobalAttributes

    ! Arguments

    character (LEN=*), INTENT(IN) :: l3StartDay, l3EndDay, tempfile
    Integer, intent(out) :: nDays
    integer, intent(in), optional :: L3MFlag
    
    ! Local variables

    type(MLSFile_T)                :: MLSFile
    character (LEN=FileNameLen) :: pcfNames(40)
    integer :: mis_Days(maxMisDays), misDInt
    integer :: numDays, mis_numDay, the_hdfVersion, file_id
    integer :: record_length, i, status, mis_numDays, n=1 
    integer :: indx, j, k, l3DayInt, m=1
    integer, external :: HE5_EHrdglatt
    character (len=20) :: name
    character (len=3) :: date, misD 
    character (len=7) :: misDchar 
    character (LEN=FileNameLen) :: pcfFileNames

    call GetL2GPfromPCF(mlspcf_l2gp_start,mlspcf_l2gp_end, tempfile, &
	& l3StartDay, l3EndDay, numDays, pcfNames, mis_numDays, mis_Days)

    if (present(L3MFlag)) then
	numDays = numDays + 1
        call Allocate_test (GlobalAttributes%OrbNumDays, max_orbits, &
	   & numDays,'GlobalAttributes%OrbNumDays', ModuleName)
        call Allocate_test (GlobalAttributes%OrbPeriodDays, max_orbits, &
	   & numDays,'GlobalAttributes%OrbPeriodDays', ModuleName)
        do i = 1, numDays-1
           the_hdfVersion = mls_hdf_version(trim(pcfNames(i)))
           status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_READ, &
            & name=trim(pcfNames(i)), HDFVersion=the_hdfVersion )
           call MLS_OpenFile( MLSFile )
           file_ID = MLSFile%FileID%f_id
           ! file_id = mls_io_gen_openF(l_swath, .TRUE., status, &
           !   & record_length, DFACC_READ, pcfNames(i), &
           !   & hdfVersion=the_hdfVersion, debugOption=.false. )
           !if ( status /= 0 ) &
           !   call MLSMessage ( MLSMSG_Error, ModuleName, &
           !     & "Unable to open L2gp file: " &
           !     & // trim(pcfNames(i))//' for reading')
              status = HE5_EHrdglatt(file_id, 'OrbitNumber', &
                & GlobalAttributes%OrbNumDays(:,i))
              if (status /= 0) then
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & "Unable to read attribute OrbitNumber for " &
                  & //trim(pcfNames(i)))
                GlobalAttributes%OrbNumDays(:,i) = -1
              end if

              status = HE5_EHrdglatt(file_id, 'OrbitPeriod', &
                & GlobalAttributes%OrbPeriodDays(:,i))
              if (status /= 0) then
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                   & "Unable to read attribute OrbitPeriod for "&
                   & //trim(pcfNames(i)))
                GlobalAttributes%OrbPeriodDays(:,i) = -1.0
              end if
         enddo

         GlobalAttributes%OrbNumDays(:,numDays) = -1
         GlobalAttributes%OrbNumDays(1:1,numDays) = &
	      &	GlobalAttributes%OrbNumDays(1:1,1)
         GlobalAttributes%OrbNumDays(2:2,numDays) = &
	      &	maxval(GlobalAttributes%OrbNumDays(:,numDays-1))
!        call output('OrbNum_last =', advance='no')
!        call output(GlobalAttributes%OrbNumDays(:,numDays), advance='yes')
         GlobalAttributes%OrbPeriodDays(:,numDays) = -1.0
         GlobalAttributes%OrbPeriodDays(1:1,numDays) = &
	      &	GlobalAttributes%OrbPeriodDays(1:1,1)
         GlobalAttributes%OrbPeriodDays(2:2,numDays) = &
	      &	maxval(GlobalAttributes%OrbPeriodDays(:,numDays-1))
!       call output('OrbNum_last =', advance='no')
!       call output(GlobalAttributes%OrbPeriodDays(:,numDays), advance='yes')
    else 
        call Allocate_test (GlobalAttributes%OrbNumDays, max_orbits, &
	   & numDays+mis_numDays,'GlobalAttributes%OrbNumDays', ModuleName)
        call Allocate_test (GlobalAttributes%OrbPeriodDays, max_orbits, &
	   & numDays+mis_numDays,'GlobalAttributes%OrbPeriodDays', ModuleName)

        DO i = 1, numDays
       	   m = n
           the_hdfVersion = mls_hdf_version(trim(pcfNames(i)))
           pcfFileNames = trim(pcfNames(i))
           indx = INDEX(pcfFileNames, '.', .TRUE.)
           date = pcfFileNames(indx-3:indx-1)
           read(date, '(i3)')l3DayInt 
           l3DayInt = l3DayInt - 1
       
           do j = 1, mis_numDays
	      write(misDchar, '(I7)') mis_Days(j)
	      misD = misDchar(5:7)
              read(misD, '(i3)')misDInt 
	      if (l3DayInt == misDInt) then
	         do k = 1, 16
	           GlobalAttributes%OrbNumDays(k,m) = &
		     & GlobalAttributes%OrbNumDays(15,m-1)+k
	           GlobalAttributes%OrbPeriodDays(k,m) = 5930.5056 
	         enddo  
	         n = n+1
	         m = m+1
	      endif
          enddo

           status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_READ, &
            & name=trim(pcfNames(i)), HDFVersion=the_hdfVersion )
           call MLS_OpenFile( MLSFile )
           file_ID = MLSFile%FileID%f_id
          ! file_id = mls_io_gen_openF('swopen', .TRUE., status, &
          ! file_id = mls_io_gen_openF(l_swath, .TRUE., status, &
          !  & record_length, DFACC_READ, pcfNames(i), &
          !  & hdfVersion=the_hdfVersion, debugOption=.false. )
          !if ( status /= 0 ) &
          !  call MLSMessage ( MLSMSG_Error, ModuleName, &
		    !& "Unable to open L2gp file: " &
		    !& // trim(pcfNames(i))//' for reading')
            status = HE5_EHrdglatt(file_id, 'OrbitNumber', &
		& GlobalAttributes%OrbNumDays(:,m))
	    n = n+1
      	    if (status /= 0) then 
               call MLSMessage ( MLSMSG_Warning, ModuleName, &
       	    	  & "Unable to read attribute OrbitNumber, calculate "  &
		  & //trim(pcfNames(i)))
	       do j = 1, 16
	          GlobalAttributes%OrbNumDays(j,m) = &
		     & GlobalAttributes%OrbNumDays(15,m-1)+j
	       enddo  
            end if

!            call output(GlobalAttributes%OrbNumDays(:,m), advance='yes')

            status = HE5_EHrdglatt(file_id, 'OrbitPeriod', & 
		& GlobalAttributes%OrbPeriodDays(:,m))
            if (status /= 0) then
         	call MLSMessage ( MLSMSG_Warning, ModuleName, &
       	    	   & "Unable to read attribute OrbitPeriod for "&
		   & //trim(pcfNames(i)))
	       do j = 1, 16
	          GlobalAttributes%OrbPeriodDays(j,m) = &
		     & GlobalAttributes%OrbPeriodDays(15,m-1)
	       enddo  
	  !	GlobalAttributes%OrbPeriodDays(:,m) = -1.0
		
            end if

!            call output(GlobalAttributes%OrbPeriodDays(:,m), advance='yes')

            status = HE5_EHrdglatt(file_id, 'FirstMAF', &
                & GlobalAttributes%FirstMAFCtr)
            if (status /= 0) then
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                   & "Unable to read attribute FirstMAF for "&
                   & //trim(pcfNames(i)))
            end if
            call output(GlobalAttributes%FirstMAFCtr, advance='yes')
                                                                            
            status = HE5_EHrdglatt(file_id, 'LastMAF', &
                & GlobalAttributes%LastMAFCtr)
            if (status /= 0) then
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                   & "Unable to read attribute LastMAF for "&
                   & //trim(pcfNames(i)))
            end if
            call output(GlobalAttributes%LastMAFCtr, advance='yes')

            call MLS_CloseFile( MLSFile )
            ! status = mls_io_gen_closeF('swclose', file_id, pcfNames(i), & 
            ! status = mls_io_gen_closeF(l_swath, file_id, pcfNames(i), & 
           	!& hdfVersion=the_hdfVersion)
            !if ( status /= 0 ) &
           	!call MLSMessage ( MLSMSG_Error, ModuleName, &
       	   !   	   & "Unable to close L2gp file:"//trim(pcfNames(i)))
        END DO
    end if

  !-------------------------------
  END SUBROUTINE ReadL2GPAttribute
  !-------------------------------

  !--------------------------------------------------------
  SUBROUTINE ReadL2GP(L2FileName, swathname, l2gpData)
  !--------------------------------------------------------
    
    ! Brief description of subroutine
    ! This routine reads an L2GP file, returning a filled data structure.
    
    ! Arguments

    CHARACTER (LEN=*), INTENT(IN) :: swathname
    CHARACTER (LEN=FileNameLen), INTENT(IN) :: L2FileName
        
    TYPE( L2GPData_T ), INTENT(OUT) :: l2gpData
   
    call ReadL2GPData(trim(L2FileName), swathname, l2gpData)
      
  !-------------------------
  END SUBROUTINE ReadL2GP
  !-------------------------
    
  !--------------------------------------------------------------------------
  SUBROUTINE SetupL2GPProd (numDays, totalIndex, l3rtemp, l3Res)
  !--------------------------------------------------------------------------
    ! Brief description of subroutine
    ! This subroutine setup l2gp Records files for a quantity over a range of days.  
    ! It returns an array of L2GPData_T structures, 
    ! with one structure per day.

    ! Arguments
      
    INTEGER, INTENT(IN) :: numDays, totalIndex
      
    TYPE( L2GPData_T ), POINTER :: l3Res(:), l3rtemp(:)
      
    ! Parameters
            
    ! Variables
      
    CHARACTER (LEN=480) :: msr

    INTEGER :: err, i, nTimes, nTimesTotal, nFreqs, nLevels

    LOGICAL :: Fillin

    IF (numDays == 0) THEN
         
       ! If not, issue a warning
         
       msr = 'No L2GP data found for '
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
       NULLIFY(l3Res)
         
       ! If so, open, read, & close the files for this quantity
         
    ELSE
         
       ALLOCATE( l3Res(numDays), STAT=err )
       IF ( err /= 0 ) THEN
          msr = MLSMSG_Allocate // ' l2gp array for '
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
        
       DO i = 1, numDays
          ! Read information from the L2GP file
          CALL SetupNewL2GPRecord(l3Res(i), nFreqs=1, nLevels=totalIndex, &
	     & nTimes=l3rtemp(i)%nTimes)
          l3Res(i)%chunkNumber = l3rtemp(i)%chunkNumber
          l3Res(i)%solarTime = l3rtemp(i)%solarTime
          l3Res(i)%solarZenith = l3rtemp(i)%solarZenith
          l3Res(i)%losAngle = l3rtemp(i)%losAngle
          l3Res(i)%geodAngle = l3rtemp(i)%geodAngle
          l3Res(i)%latitude = l3rtemp(i)%latitude
          l3Res(i)%longitude = l3rtemp(i)%longitude
          l3Res(i)%time = l3rtemp(i)%time
          l3Res(i)%quality = l3rtemp(i)%quality
          l3Res(i)%status = l3rtemp(i)%status
          l3Res(i)%frequency = l3rtemp(i)%frequency
          l3Res(i)%nTimes = l3rtemp(i)%nTimes
          l3Res(i)%nTimesTotal = l3rtemp(i)%nTimesTotal
          l3Res(i)%nLevels = totalIndex
          l3Res(i)%nFreqs = 0
          l3Res(i)%verticalCoordinate = 'Pressures'
       ENDDO
            
    ENDIF
      
    !-----------------------------
    END SUBROUTINE SetupL2GPProd
    !-----------------------------
  !--------------------------------------------------------------------------
  SUBROUTINE ReadL2GPProd (product, template, startDOY, endDOY, numDays, & 
       & mis_numDays, mis_Days, prodL2GP)
  !--------------------------------------------------------------------------
      
    ! Brief description of subroutine
    ! This subroutine reads l2gp files for a quantity over a range of days.  
    ! It returns an array of L2GPData_T structures, 
    ! with one structure per day.

    ! Arguments
      
    CHARACTER (LEN=*), INTENT(IN) :: product, template
      
    CHARACTER (LEN=8), INTENT(IN) :: endDOY, startDOY
      
    INTEGER, INTENT(OUT) :: numDays, mis_numDays
      
    TYPE( L2GPData_T ), POINTER :: prodL2GP(:)
      
    !CHARACTER (LEN=DATE_LEN), INTENT(OUT) :: mis_Days(:)
    INTEGER, INTENT(OUT) :: mis_Days(:)
      
    ! Parameters
            
    ! Variables
    CHARACTER (LEN=480) :: msr
    CHARACTER (LEN=FileNameLen) :: pcfNames(40)
    CHARACTER (LEN=8) :: date

    INTEGER :: err, i
    
    mis_Days = 0
    mis_numDays = 0        
    ! Search the PCF for L2GP files for this species
      
    CALL GetL2GPfromPCF(mlspcf_l2gp_start, mlspcf_l2gp_end, template, & 
         & startDOY, endDOY, numDays, pcfNames, mis_numDays, mis_Days)
      
    ! Check that numDays is at least 1, so pointer allocation is possible
      
    IF (numDays == 0) THEN
         
       ! If not, issue a warning
         
       msr = 'No L2GP data found for ' // product
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
       NULLIFY(prodL2GP)
         
       ! If so, open, read, & close the files for this quantity
         
    ELSE
         
       ALLOCATE( prodL2GP(numDays), STAT=err )
       IF ( err /= 0 ) THEN
          msr = MLSMSG_Allocate // ' l2gp array for ' // product
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       DO i = 1, numDays
          ! Read information from the L2GP file
          CALL ReadL2GP(pcfNames(i), product, prodL2GP(i) )
       ENDDO
            
    ENDIF
      
    !-----------------------------
    END SUBROUTINE ReadL2GPProd
    !-----------------------------
 
 !---------------------------------------
 SUBROUTINE ResidualOutput (type, l3r)
 !---------------------------------------
  USE HDF, ONLY: DFACC_RDWR
   
   ! Brief description of subroutine
   ! This program creates & writes the L3Residual swaths in an l3 map file.
   
   ! Arguments
   
   CHARACTER (LEN=*), INTENT(IN) :: type
   
   TYPE( L2GPData_T ), INTENT(INOUT) :: l3r(:)

   ! Functions
         
   ! Parameters
   
   ! Variables
   
   type(MLSFile_T)                :: MLSFile
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: l3File
   CHARACTER (LEN=8) :: date
   REAL(r8) :: validTime
   INTEGER :: i, j, match, numDays, file_id, hdfVersion, record_length, status 
   LOGICAL :: notUnlimited 

   ! For each element of the database,
   
   numDays = SIZE(l3r)
   
   DO i = 1, numDays
      
      ! Find the l3dm file for the proper day
   
      DO j = 1, l3r(i)%ntimes
         validTime = l3r(i)%time(j)
         IF (validTime /= 0.0) Exit
      ENDDO
  
      CALL FindFileDay(type, validTime, mlspcf_l3dm_start, &
           & mlspcf_l3dm_end, match, l3File, date)
      IF (match == -1) THEN
         msr = 'No ' // TRIM(type) // ' file for day ' // date & 
              & // ' found for appending ' // TRIM(l3r(i)%name) // 'swath.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      
      ! Re-open the l3dm grid file for the creation of swaths
      ! Create the L3Residual swath in this day's file
      ! Re-attach to the swath for writing
      ! Write to the geolocation & data fields

      hdfVersion = mls_hdf_version(trim(l3File))

      status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_RDWR, &
       & name=l3File, HDFVersion=hdfVersion )
      call MLS_OpenFile( MLSFile )
      file_ID = MLSFile%FileID%f_id
      ! file_id = mls_io_gen_openF('swopen', .TRUE., status, &
      ! file_id = mls_io_gen_openF(l_swath, .TRUE., status, &
      !     & record_length, DFACC_RDWR, l3File, &
      !     & hdfVersion=hdfVersion, debugOption=.false. )

      ! call dump (real(l3r(i)%l2gpValue, r8), 'l2gpValue: ')
      call WriteL2GPData(l3r(i), file_id, l3r(i)%name, & 
           & hdfVersion=hdfVersion,notUnlimited=.true.)

      ! status = mls_io_gen_closeF(l_swath, file_id, l3File, & 
      !     & hdfVersion=hdfVersion)
      call MLS_CloseFile ( MLSFile )

   ENDDO
      
 !-------------------------------
 END SUBROUTINE ResidualOutput
 !-------------------------------

 !-----------------------------------------------------------------------
 SUBROUTINE ReadL2DGData (product, numFiles, files, numDays, prodL2DG)
   !-----------------------------------------------------------------------
      
   ! Brief description of subroutine
   ! This subroutine reads l2gp dg files for a quantity over range of days.
   ! It returns an array of L2GPData_T structures, 
   ! with one structure per day.
   
   ! Arguments
      
   CHARACTER (LEN=*), INTENT(IN) :: product
   
   CHARACTER (LEN=FileNameLen), INTENT(IN) :: files(:)
      
   INTEGER, INTENT(IN) :: numFiles
   
   INTEGER, INTENT(OUT) :: numDays
      
   TYPE( L2GPData_T ), POINTER :: prodL2DG(:)

   ! Functions
   INTEGER, EXTERNAL :: he5_swattach, he5_swclose, &
           & he5_swdetach, he5_swopen, he5_swinqswath
                                                                            
   ! Parameters
   
   ! Variables
      
   CHARACTER (LEN=12000) :: list, templist
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: prodFiles(numFiles)
   CHARACTER (LEN=GridNameLen) :: swath
      
   INTEGER :: err, i, j, indx, len, ns, hdfversion
   
   ! Check for the product in the input files
      
   numDays = 0
  
   ns = he5_swinqswath( trim(files(1)), templist, len)
                                                                            
   IF (ns == -1) THEN
         msr = 'Failed to read swath list from ' // files(1)
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         !CYCLE
   ENDIF

   DO i = 1, numFiles
      list = templist 
      DO j = 1, ns
            
         indx = INDEX(list, ',')
         
         IF (indx == 0) THEN
            swath = list
         ELSE
            swath = list(:indx-1)
            list  = list(indx+1:)

         ENDIF
            
         IF ( TRIM(swath) == TRIM(product) ) THEN
            numDays = numDays + 1
            prodFiles(numDays) = files(i)
            EXIT
         ENDIF
      ENDDO
   ENDDO
      
   ! Check that numDays is at least 1, so pointer allocation is possible

   IF (numDays == 0) THEN
         
      ! If not, issue a warning
      
      msr = 'No L2GP data found for ' // product
      CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
      NULLIFY(prodL2DG)
         
      ! If so, open, read, & close the files for this quantity
         
   ELSE

      ALLOCATE( prodL2DG(numDays), STAT=err )
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' l2gp array for ' // product
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
         
      DO i = 1, numDays
         
         ! Read information from the L2GP file
         
         CALL ReadL2GP( prodFiles(i), product, prodL2DG(i) )
                        
      ENDDO
         
   ENDIF

   !-----------------------------
 END SUBROUTINE ReadL2DGData
 !-----------------------------

 !=====================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
END MODULE L2Interface
!=====================

!# $Log$
!# Revision 1.24  2008/12/02 23:13:56  pwagner
!# mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!#
!# Revision 1.23  2007/06/26 19:06:26  cvuu
!# fix bug
!#
!# Revision 1.22  2006/05/31 21:05:35  cvuu
!# Check missing dates that are not in the PCF; Fix orbit number and orbit period for L3DZ
!#
!# Revision 1.21  2006/05/24 19:14:30  cvuu
!# Fix the deallocate problem
!#
!# Revision 1.20  2006/04/17 15:43:04  cvuu
!# Remove the print statements
!#
!# Revision 1.19  2006/02/28 17:56:56  cvuu
!# V2.00 commit
!#
!# Revision 1.18  2005/06/23 19:07:38  pwagner
!# Reworded Copyright statement, moved rcs id
!#
!# Revision 1.17  2005/06/14 20:46:42  pwagner
!# Changed interface to mls_io_gen functions
!#
!# Revision 1.16  2004/06/22 16:55:47  cvuu
!# Using SetupNewL2GPRecord for L3 Residual based on the nLevels from L3 not L2
!#
!# Revision 1.15  2004/05/04 15:33:15  cvuu
!# v1.4.3: Use int array for Date in Data Field
!#
!# Revision 1.14  2004/01/07 21:43:18  cvuu
!# version 1.4 commit
!#
!# Revision 1.13  2003/09/15 18:30:48  cvuu
!# Read OrbitNumber and OrbitPeriod from L2GP files
!#
!# Revision 1.12  2003/04/30 18:15:48  pwagner
!# Work-around for LF95 infinite compile-time bug
!#
!# Revision 1.11  2003/03/22 02:41:32  jdone
!# implemented L2GPData library routines
!#
!# Revision 1.10  2002/04/15 22:25:00  jdone
!#  swid checked
!#
!# Revision 1.9  2002/04/11 22:29:47  jdone
!# swid checked
!#
!# Revision 1.8  2002/04/10 21:54:05  jdone
!# error checking on allocated statements
!#
!# Revision 1.7  2002/02/20 19:21:15  ybj
!# *** empty log message ***
!#
!# Revision 1.6  2001/07/18 15:50:59  nakamura
!# Generalized to work with L3M as well.
!#
!# Revision 1.5  2001/04/24 19:37:30  nakamura
!# Changes for privatization of L2GPData.
!#
!# Revision 1.4  2001/02/21 20:37:18  nakamura
!# Changed MLSPCF to MLSPCF3.
!#
!# Revision 1.3  2000/12/29 20:40:04  nakamura
!# Changed ReadL2GPProd to take start & end days as input; modified ResidualOutput for the one-product/all-days paradigm.
!#
!# Revision 1.2  2000/12/07 19:28:27  nakamura
!# Updated for level becoming an expandable template field.
!#
!# Revision 1.1  2000/11/15 20:56:51  nakamura
!# Module containing subroutines related to L2GP data.
!#
