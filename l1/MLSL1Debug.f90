! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSL1Debug! Module for debugging L1 code
!=============================================================================

!+
!
! A module to use debugging MLS Level 1 code. Can be used to dump data while
! processing.
!
! To use set the env variable MLSL1DEBUG.
!
! sh-like shells
!
! % export MLSL1DEBUG=list,of,filenames,that,turn,on,debugging
!
! csh-like shells
!
! % setenv MLSL1DEBUG list,of,filenames,that,turn,on,debugging
!
! Each element of the list is the name of a file to write debugging information
! to. The beginning of the filename tells which debugging code to turn
! on. The filename can be anything after the initial substring -- which is what
! this module keys on -- but the first substring is required.
!
! Currently, you have the following choices. (the comparisons are case
! insensitive, but the string you pass is used *as is*)
!
! 'comvec...' : writes information about the comvecs computation from
!               Calibration.f90. INitializes in InitCalibWindow and writes info
!               out of SetComVecs. The filanem
!
! 'radiances...' : 
!
!
!-

  USE MLSL1RunConfig, ONLY :  MLSL1Executable ! name of executable currently
                                              ! running. Set in each main

  USE IO_STUFF, ONLY: get_lun

  USE MLSL1Common, ONLY: R8,FBchans, FBnum, MBchans, MBnum, &
       WFchans, WFnum, DACSchans, DACSnum, &
       GHzNum, MaxAlts, MaxMIFs, deflt_gain, BandWidth, L1ProgType,THzType


  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

  USE MLSStrings, ONLY : streq
  USE MLSStringLists, ONLY: SwitchDetail, NumStringElements,List2Array

  IMPLICIT NONE
  

  PRIVATE 

  !! Type that controls the `file' type debugging. If MLSL1DEBUG has the strings
  !! 'comvecs', 'radiances' or 'interps' *and* none of those have 'snoop' in
  !! them, then the user is requesting we write some information out. Which we
  !! write out (and it can be any or all of them) depends on what's in
  !! MLSL1DEBUG. If MLSL1DEBUG='comvecs.date,radiances.dat,interpfile.dat' it
  !! will write out all three datasets.
  TYPE MLSL1DebugCntrlLogicals_T
    LOGICAL :: ComVecs=.FALSE.  ! true when MLSL1DEBUG has 'comvecs' 
    LOGICAL :: Radiances=.FALSE. !true when MLSL1DEBUG has 'radiances' 
    LOGICAL :: Interps=.FALSE. ! true when MLSL1DEBUG has 'interps'
  END TYPE MLSL1DebugCntrlLogicals_T


  TYPE MLSL1DebugFileInfo_T
     CHARACTER(len=512) :: filename
     INTEGER(4):: lun
     LOGICAL :: IsOpen=.FALSE.
  END TYPE MLSL1DebugFileInfo_T


  TYPE (MLSL1DebugCntrlLogicals_T), SAVE ::  DebugControl
  TYPE (MLSL1DebugFileInfo_T), SAVE :: ComVecsFileInfo, RadiancesFileInfo, InterpsFileInfo


  LOGICAL :: initialized=.FALSE., SnoopComVec=.FALSE., SnoopRadiances=.FALSE.

  INTEGER :: altIndex, dbMAF, dbStartIndex,dbEndIndex
  REAL(r8) dbTAI93
  CHARACTER(len=2)comVecRecType ! either 'LS','S', or 'T'. 


  Integer:: maxCalSize ! same as Calibration::max_cal_index, but I can't get to
                       ! it because of circularity, so I have to use _Init to
                       ! set this.

  PUBLIC :: DebugControl, dbMAF, dbStartIndex,dbEndIndex,dbTAI93, comVecRecType,&
       &    SnoopRadiances, SnoopComVec
  PUBLIC :: MLSL1Debug_Init,MLSL1DebugDoDebug, &
       &    openMLSL1DebugFiles,closeMLSL1DebugFiles, &
       &    comVecInfo_init,writeComVecInfo, writeRadiancesInfo, writeInterpsInfo



!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  CONTAINS

    ! Initialize the module by checking to see if this run has requested that we
    ! execute the debug code.

    ! ==================================================================================================
    SUBROUTINE MLSL1Debug_Init()

      ! Initialize the module by setting the debug and initialized logicals
      DebugControl%ComVecs=.FALSE.
      DebugControl%Radiances=.FALSE.
      DebugControl%Interps=.FALSE.

      ! set file info structures
      RadiancesFileInfo%lun=-1
      RadiancesFileInfo%isOpen=.FALSE.

      ComVecsFileInfo%isOpen=.FALSE.
      ComVecsFileInfo%lun=-1

      InterpsFileInfo%lun=-1
      InterpsFileInfo%isOpen=.FALSE.

      
      CALL MLSL1DebugDoDebug()
      initialized=.TRUE.
    END SUBROUTINE mlsl1debug_init

    ! ==================================================================================================
    SUBROUTINE MLSL1DebugDoDebug()
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxEls=10
      CHARACTER(len=256), dimension(maxEls):: list
      CHARACTER(len=512) :: msg 
      INTEGER :: strLen

      CHARACTER(len=512) :: value
      INTEGER:: status
      INTEGER:: ii, nn ! counters

      ! find if the MLSL1DEBUG environmental variable is set. This will be a
      ! list of one element or a comma separated list. The elements of the list
      ! are the names of the files to which the data for each debuggable module
      ! will be written, or it can be a list of words beginning with `snoop' to
      ! tell the executable which modules to `snoop'. Currently implemented is
      ! `snoopradiances'
      ! 

      IF (L1ProgType .EQ. THzType) THEN 
         PRINT *,"We don't do THz at the moment: no Debugging available!"
         return
      ENDIF
      strLen=512
      CALL GET_ENVIRONMENT_VARIABLE('MLSL1DEBUG',value,strLen,status,.TRUE.)
      IF (status==1) THEN
        ! variable isn't defined, so clearly we're not going to be doing any
        ! debugging
        PRINT *,"Env Var `MLSL1DEBUG' isn't defined! No debugging will be done!"
        CALL MLSMessage (MLSMSG_Info, ModuleName, &
             & "Env Var MLSL1DEBUG isn't defined: No debugging!")

        return
      ELSE IF (status/=0) THEN 
        PRINT *,"Failure getting MLSL1DEBUG env var"
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
             & "Failure getting MLSL1DEBUG env var")
      ENDIF

      IF (TRIM(VALUE)=="") THEN 
        CALL MLSMessage (MLSMSG_Info, ModuleName, &
             & "Env Var is empty, No debugging!")
        print *,"Env Var is empty, No debugging!"
        RETURN  ! empty env variable, nothing to do
      ENDIF

      WRITE(msg,'(a,a)') 'MLSL1DEBUG = ',trim(VALUE)
      PRINT *,TRIM(msg)
      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & TRIM(msg))

      nn=NumStringElements(TRIM(VALUE),.TRUE.,',')
      IF (nn>maxEls) THEN
        WRITE(msg,'(a,1x,i3)') &
             & 'Too many elements in MLSL1DEBUG, max allowed = ',maxEls
        PRINT *,TRIM(msg)
        CALL MLSMessage (MLSMSG_Error, ModuleName, &
             & TRIM(msg))
      ENDIF

      ! split into an array of values
      CALL List2Array(TRIM(VALUE), list, .FALSE.,',',.TRUE.)

      ! ComVec debugging
      IF (strEQ(TRIM(value),'*comvec*','wcf')) THEN
        IF (strEQ(TRIM(VALUE),'*snoop*','wcf')) THEN
          SnoopComvec=.TRUE.
          WRITE(msg,'(a)') 'Snooping Comvec!'
          CALL MLSMessage (MLSMSG_Warning, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)

        ELSE 
          ! User wants to write data to a file. Find the filename to use.
          DO ii=1,nn
            IF (strEQ(list(ii),'comvec*','wcf') ) THEN 
              DebugControl%ComVecs=.TRUE.
              ComVecsFileInfo%filename=TRIM(list(ii))
              EXIT
            ENDIF
          END DO
        ENDIF
      ENDIF

      ! Radiances debugging
      IF (strEQ(TRIM(VALUE),'*rad*','wcf') ) THEN 
        IF ( strEQ(TRIM(VALUE),'*snoop*','wcf')) THEN
          SnoopRadiances=.TRUE.
          WRITE(msg,'(a)') 'Snooping Radiances!'
          CALL MLSMessage (MLSMSG_Warning, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)

        ELSE
          ! User wants to write Radiance info to a file. Find the filename to use
          DO ii=1,nn
            IF (strEQ(list(ii),'rad*','wcf') ) THEN 
              DebugControl%Radiances=.TRUE.
              RadiancesFileInfo%filename=TRIM(list(ii))
              EXIT
            ENDIF
          END DO
        ENDIF
      ENDIF

      ! interp debugging
      IF (strEQ(TRIM(VALUE),'*interp*','wcf') ) THEN 
        IF ( strEQ(TRIM(VALUE),'*snoop*','wcf')) THEN
          SnoopRadiances=.TRUE.
          IF (snoopRadiances) THEN
            WRITE(msg,'(a)') 'Snooping interps is not implemented yet!'
            WRITE(msg,'(a)') 'Turning it off!'
          SnoopRadiances=.FALSE.
          ENDIF
          ! WRITE(msg,'(a)') 'Snooping interps!'
          ! CALL MLSMessage (MLSMSG_Warning, ModuleName, &
          !      & TRIM(msg))
          ! PRINT *,TRIM(msg)

        ELSE
          ! User wants to write the INTERP info to a file. Find the filename to use
          DO ii=1,nn
            IF (strEQ(list(ii),'interp*','wcf') ) THEN 
              DebugControl%Interps=.TRUE.
              InterpsFileInfo%filename=TRIM(list(ii))
              EXIT
            ENDIF
          END DO
        ENDIF
      ENDIF
    END SUBROUTINE MLSL1DebugDoDebug

    ! ==================================================================================================
    
    SUBROUTINE openMLSL1DebugFiles()

      ! Open one or all of the files currently supported by this module
      IMPLICIT NONE 
      CHARACTER(len=512)  ::  msg
      logical:: exists
      INTEGER:: lun
      INTEGER:: iostat=0
      CHARACTER(len=512) :: filename

      ! does the user want to debug ComVec calculations? (Calibration.f90::SetComVecs)
      IF (DebugControl%ComVecs) THEN
        filename=ComVecsFileInfo%filename
        call MLSMessage(MLSMSG_Info,ModuleName,&
             & "Opening "//TRIM(filename))
        print *,"Opening "//TRIM(filename)
        INQUIRE(file=TRIM(filename),exist=exists)
        IF (exists) THEN 
          WRITE(msg,'(a,a)') TRIM(filename),' exists! Overwriting!'
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        
        CALL get_lun(lun)
        IF (lun == -1) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't get LUN to open ",TRIM(filename),",lun=",&
               & lun
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF

        OPEN(unit=lun,file=TRIM(filename),&
             & access='stream',form='unformatted',&
             & status='replace', ioStat=iostat)
        IF (iostat .NE. 0) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't open ",TRIM(filename),",iostat=",&
               & iostat
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        ComVecsFileInfo%lun=lun
        ComVecsFileInfo%IsOpen=.TRUE.
      ENDIF ! debug ComVecs stuff


      ! does the user want to debug radiances? (Radiances::CalcRadiances)
      IF (DebugControl%Radiances) THEN
        filename=RadiancesFileInfo%filename
        call MLSMessage(MLSMSG_Info,ModuleName,&
             & "Opening "//TRIM(filename))
        print *,"Opening "//TRIM(filename)
        INQUIRE(file=TRIM(filename),exist=exists)
        IF (exists) THEN 
          WRITE(msg,'(a,a)') TRIM(filename),' exists! Overwriting!'
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        
        CALL get_lun(lun)
        IF (lun == -1) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't get LUN to open ",TRIM(filename),",lun=",&
               & lun
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF

        OPEN(unit=lun,file=TRIM(filename),&
             & access='stream',form='unformatted',&
             & status='replace', ioStat=iostat)
        
        IF (iostat .NE. 0) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't open ",TRIM(filename),",iostat=",&
               & iostat
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        RadiancesFileInfo%lun=lun
        RadiancesFileInfo%IsOpen=.TRUE.
      ENDIF ! test on debug Radiances stuff

      IF (DebugControl%Interps) THEN
        filename=InterpsFileInfo%filename
        call MLSMessage(MLSMSG_Info,ModuleName,&
             & "Opening "//TRIM(filename))
        print *,"Opening "//TRIM(filename)
        INQUIRE(file=TRIM(filename),exist=exists)
        IF (exists) THEN 
          WRITE(msg,'(a,a)') TRIM(filename),' exists! Overwriting!'
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        
        CALL get_lun(lun)
        IF (lun == -1) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't get LUN to open ",TRIM(filename),",lun=",&
               & lun
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF

        OPEN(unit=lun,file=TRIM(filename),&
             & access='stream',form='unformatted',&
             & status='replace', ioStat=iostat)
        IF (iostat .NE. 0) THEN 
          WRITE(msg,'(a,a,a,i4)') &
               & "Couldn't open ",TRIM(filename),",iostat=",&
               & iostat
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & TRIM(msg))
          PRINT *,TRIM(msg)
        ENDIF
        InterpsFileInfo%lun=lun
        InterpsFileInfo%IsOpen=.TRUE.
      ENDIF ! debug Interps stuff

    END SUBROUTINE openMLSL1DebugFiles

    ! ==================================================================================================
    ! close all files
    SUBROUTINE closeMLSL1DebugFiles()

      IF (DebugControl%ComVecs .AND. ComVecsFileInfo%isOpen) THEN
        CALL MLSMessage(MLSMSG_Info,ModuleName,"Closing ComVecsFile")
        print *,"Closing ComVecsFile"
        CLOSE(ComVecsFileInfo%lun)
        ComVecsFileInfo%IsOpen=.FALSE.
      ENDIF ! test on debug ComVecs stuff

      IF (DebugControl%Radiances .AND. RadiancesFileInfo%isOpen) THEN
        CALL MLSMessage(MLSMSG_Info,ModuleName,"Closing RadiancesFile")
        print *,"Closing RadiancesFile"

        CLOSE(RadiancesFileInfo%lun)
        RadiancesFileInfo%IsOpen=.FALSE.
      ENDIF ! test on debug Radiances stuff

      IF (DebugControl%Interps .AND. InterpsFileInfo%isOpen) THEN
        CALL MLSMessage(MLSMSG_Info,ModuleName,"Closing InterpsFile")
        print *,"Closing InterpsFile"
        CLOSE(InterpsFileInfo%lun)
        InterpsFileInfo%IsOpen=.FALSE.
      ENDIF ! test on debug Interps stuff

    END SUBROUTINE closeMLSL1DebugFiles

    ! ==================================================================================================
    ! Initialize the code used for debugging Calibration::ComVec processing


    SUBROUTINE comVecInfo_Init(retType,maxMIFS,maxAlts,max_cal_index)
      !! Writes a header record to the ComVec debug file giving run information,
      !! mostly having to do with sizing info, so that the reader can then use
      !! that info to dimension arrays correctly. WinMAFs, MaxMIFs,
      !! max_cal_index

      CHARACTER(len=2),intent(in) :: retType
      INTEGER,intent(in):: maxMIFS,maxAlts,max_cal_index
      maxCalSize=max_cal_index

      write(ComVecsFileInfo%Lun) comVecRecType, maxMIFs,maxAlts,max_cal_index
      
    END SUBROUTINE comVecInfo_init

    ! ==================================================================================================
    ! Write a record for ComVecInfo processing
    SUBROUTINE writeComVecInfo(altIndex,MAF,TAI93,seq,qual,weight,time,errmul,comvec)

      INTEGER, intent(in)          :: altIndex
      REAL(r8),intent(in)          :: TAI93
      CHARACTER(len=1), intent(in) :: seq(:)
      INTEGER,INTENT(in)           :: MAF, qual(:) 
      REAL(r8),intent(in)          :: weight(:)  
      REAL(r8),intent(in)          :: time(:)    
      REAL(r8),INTENT(in)          :: errmul(:)
      REAL(r8),INTENT(in)          :: comvec(:,:)

      IF (DebugControl%ComVecs .EQV. .FALSE.) RETURN

      WRITE(ComVecsFileInfo%Lun) comVecRecType,  &
           &                    altIndex,&   
           &                    MAF, &
           &                    tai93, &     
           &                    seq,&        
           &                    qual,&       
           &                    weight,&     
           &                    time,&       
           &                    errmul,&     
           &                    comvec       
      ! CALL FLUSH(ComVecsFileInfo%Lun)
    END SUBROUTINE writeComVecInfo

    ! =========================================================================================
    ! Setup writing Radiances
    ! =========================================================================================
    
    SUBROUTINE radiancesInfo_Init(max_cal_index)
      INTEGER,intent(in):: max_cal_index
      maxCalSize=max_cal_index
    END SUBROUTINE  radiancesInfo_Init


    ! =========================================================================================
    ! Write Radiances information to the radiances debug file
    ! =========================================================================================
    SUBROUTINE writeRadiancesInfo(gnum,nchan,MAFNo, counterMAF,&
         & sind,eind, center, target_t, space_t, GHz_T1,GHz_T2)


      USE MLSL1Common, ONLY: limb_cnts, space_interp, target_interp, &
           &                   target_err, slimb_interp, slimb_err, &
           &                   slimb_type, space_err

      USE MLSL1Rad, ONLY : FBrad

      USE BandTbls, ONLY : SpilloverLoss


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAFno, counterMAF,gnum, nchan, sind, eind, center
      REAL(4), INTENT(IN) :: target_t, space_t, GHz_T1,GHz_T2


      ! Local variables 
      INTEGER :: eind2,i,c,b ! counters

      eind2=eind-sind
      IF (DebugControl%Radiances .EQV. .FALSE.) RETURN
      print *,'Writing info to the Radiances file'
      WRITE(RadiancesFileInfo%Lun) CounterMAF,MAFno, center, &
           &                       gnum,nchan,sind,eind,target_t, space_t,&
                     &             GHz_T1,GHz_T2,                                               &
                     &             (FBRad(i)%bandno, i=1,gnum),                                 &
                     &             ((SpilloverLoss(b)%eta_TSL(c,1),b=1,gnum),c=1,3),            &
                     &             (FBRad(i)%signal%radiometerNumber, i=1,gnum),                &
                     &            (((limb_cnts%FB(i,c,b), i=sind,eind),c=1,nchan),b=1,gnum),    &
                     &            ((BandWidth%FB(c,b) , c=1,nchan),b=1,gnum),                   &
                     &            (((space_interp(i)%FB(c,b),  i=0,eind2),c=1,nchan),b=1,gnum), &
                     &            ((deflt_gain%FB(c,b), c=1,nchan),b=1,gnum),                   &
                     &            (((target_interp(i)%FB(c,b), i=0,eind2),c=1,nchan),b=1,gnum), &
                     &            (((space_err(i)%FB(c,b), i=0,eind2),c=1,nchan),b=1,gnum),     &
                     &            ((slimb_type%FB(c,b), c=1,nchan),b=1,gnum),                   &
                     &            (((target_err(i)%FB(c,b), i=0,eind2),c=1,nchan),b=1,gnum),    &
                     &            (((slimb_interp(i)%FB(c,b),i=0,eind2),c=1,nchan),b=1,gnum),   &
                     &             (FBRad(i)%signal%radiometerModifier, i=1,gnum)


      
      ! CALL FLUSH(RadiancesFileInfo%Lun)
    END SUBROUTINE writeRadiancesInfo


    SUBROUTINE writeInterpsInfo(numMIFs, &
         &                      tai93,&           
         &                      space_interp,&    
         &                      space_err,&       
         &                      slimb_interp,&    
         &                      slimb_err,&       
         &                      target_interp,&   
         &                      target_err)       

      INTEGER, INTENT(in) :: numMIFs
      REAL(r8), INTENT(in) :: tai93
      REAL(r8), INTENT(in) :: space_interp(0:numMIFs),&
           &                  space_err(0:numMIFs),&
           &                  slimb_interp(0:numMIFs),& 
           &                  slimb_err(0:numMIFs),&    
           &                  target_interp(0:numMIFs),&
           &                  target_err(0:numMIFs)     

      WRITE(InterpsFileInfo%Lun) numMIFs
      WRITE(InterpsFileInfo%Lun) tai93
      WRITE(InterpsFileInfo%Lun) space_interp
      WRITE(InterpsFileInfo%Lun) space_err
      WRITE(InterpsFileInfo%Lun) slimb_interp 
      WRITE(InterpsFileInfo%Lun) slimb_err    
      WRITE(InterpsFileInfo%Lun) target_interp
      WRITE(InterpsFileInfo%Lun) target_err     

    END SUBROUTINE writeInterpsInfo      

    
END MODULE MLSL1Debug


