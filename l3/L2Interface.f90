
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!==============================================================================
MODULE L2Interface
!==============================================================================

  USE L2GPData, ONLY: L2GPData_T, ReadL2GPData, WriteL2GPData
  USE MLSCommon, ONLY: r8
  USE MLSFiles, ONLY: mls_openFile, mls_closeFile, & 
       & mls_hdf_version, mls_inqswath, mls_io_gen_openF, mls_io_gen_closeF
  USE MLSL3Common, ONLY: DATE_LEN, CCSDS_LEN, FILENAMELEN, maxWindow, & 
       & GRIDNAMELEN
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Allocate, &
       & MLSMSG_WARNING, MLSMSG_INFO, MLSMSG_FILEOPEN
  USE MLSPCF3, ONLY: mlspcf_l2gp_start, mlspcf_l2gp_end, &
       & mlspcf_l3dm_start, mlspcf_l3dm_end
  USE PCFModule, ONLY: FindFileDay, ExpandFileTemplate
  IMPLICIT NONE
  private
  PUBLIC :: GetL2GPfromPCF, ReadL2GP, ReadL2GPProd, ResidualOutput, ReadL2DGData
  
  PRIVATE :: ID, ModuleName
  
  !------------------- RCS Ident Info -----------------------
  CHARACTER(LEN=130) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !----------------------------------------------------------
  
  ! Contents:
  
  ! Subroutines -- GetL2GPfromPCF
  !                ReadL2GP
  !                ReadL2GPProd
  !                ResidualOutput
  !                ReadL2DGData
  
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
    
    CHARACTER (LEN=FileNameLen) :: pcfNames(:)
    CHARACTER (LEN=DATE_LEN), INTENT(OUT) :: mis_Days(:)
    
    ! Parameters
    
    ! Functions
    
    ! Variables
    
    CHARACTER (LEN=FileNameLen) :: physicalFilename, type
    CHARACTER (LEN=8) :: date
    
    INTEGER :: i, indx, returnStatus, version
    
    ! Expand the level/species template
    
    CALL ExpandFileTemplate(template, type, 'L2GP')
    
    numFiles = 0
    mis_numFiles = 0
    
    ! Loop through all the PCF numbers for L2GP files
    
    DO i = mlspcf_start, mlspcf_end
       
       version = 1
       returnStatus = Pgs_pc_getReference(i, version, physicalFilename)
       
       ! If no file name was returned, go on to the next PCF number
       
       IF(returnStatus /= PGS_S_SUCCESS) THEN
          
          IF ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN
             
             ! Extract the date from the file name
             
             indx = INDEX(physicalFilename, '.', .TRUE.)
             date = physicalFilename(indx-8:indx-1)
             
             ! Check that the date is within the desired boundaries for reading
             
             IF ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN
                
                ! Save the name in an array of files for the species
                
                mis_numFiles = mis_numFiles + 1
                mis_Days(mis_numFiles) = date
                
             ENDIF
             
          ENDIF
          
          CYCLE
          
       ENDIF
       
       ! Check that the returned file name is for the proper species
       
       IF ( INDEX(physicalFilename, TRIM(type)) /= 0 ) THEN
          
          ! Extract the date from the file name
          
          indx = INDEX(physicalFilename, '.', .TRUE.)
          date = physicalFilename(indx-8:indx-1)
          
          ! Check that the date is within the desired boundaries for reading
          
          IF ( (date .GE. startDOY) .AND. (date .LE. endDOY) ) THEN
             
             ! Save the name in an array of files for the species
             
             numFiles = numFiles + 1
             pcfNames(numFiles) = physicalFilename
             
          ENDIF
          
       ENDIF
       
    ENDDO
    
  !-------------------------------
  END SUBROUTINE GetL2GPfromPCF
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
    
    call ReadL2GPData(L2FileName, swathname, l2gpData)
      
  !-------------------------
  END SUBROUTINE ReadL2GP
  !-------------------------
    
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
      
    CHARACTER (LEN=DATE_LEN), INTENT(OUT) :: mis_Days(maxWindow)
      
    ! Parameters
            
    ! Variables
      
    CHARACTER (LEN=480) :: msr
    CHARACTER (LEN=FileNameLen) :: pcfNames(40)
      
    INTEGER :: err, i
            
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
   
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: l3File
   CHARACTER (LEN=8) :: date
   
   INTEGER :: i, match, numDays, file_id, hdfVersion, record_length, status 

   ! For each element of the database,
   
   numDays = SIZE(l3r)
   
   DO i = 1, numDays
      
      ! Find the l3dm file for the proper day
      
      CALL FindFileDay(type, l3r(i)%time(1), mlspcf_l3dm_start, &
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

      file_id = mls_io_gen_openF('swopen', .TRUE., status, &
           & record_length, DFACC_RDWR, l3File, &
           & hdfVersion=hdfVersion, debugOption=.false. )

      call WriteL2GPData(l3r(i), file_id, l3r(i)%name, & 
           & hdfVersion=hdfVersion)

      status = mls_io_gen_closeF('swclose', file_id, l3File, & 
           & hdfVersion=hdfVersion)

   ENDDO
      
 !-------------------------------
 END SUBROUTINE ResidualOutput
 !-------------------------------

 !-----------------------------------------------------------------------
 SUBROUTINE ReadL2DGData (product, numFiles, files, numDays, prodL2GP)
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
      
   TYPE( L2GPData_T ), POINTER :: prodL2GP(:)

   ! Functions
            
   ! Parameters
   
   ! Variables
      
   CHARACTER (LEN=6000) :: list
   CHARACTER (LEN=480) :: msr
   CHARACTER (LEN=FileNameLen) :: prodFiles(numFiles)
   CHARACTER (LEN=GridNameLen) :: swath
      
   INTEGER :: err, i, j, indx, len, ns

   
   ! Check for the product in the input files
      
   numDays = 0
   
   DO i = 1, numFiles
         
      ns = mls_inqswath( files(i), list, len )

      IF (ns == -1) THEN
         msr = 'Failed to read swath list from ' // files(i)
         CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
         CYCLE
      ENDIF
         
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
         
         CALL ReadL2GP( prodFiles(i), product, prodL2GP(i) )
                        
      ENDDO
         
   ENDIF

   !-----------------------------
 END SUBROUTINE ReadL2DGData
 !-----------------------------

 !=====================
END MODULE L2Interface
!=====================

!# $Log$
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
