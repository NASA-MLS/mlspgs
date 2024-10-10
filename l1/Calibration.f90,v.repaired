! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE Calibration ! Calibration data and routines
!=============================================================================

  USE MLSL1Debug, ONLY: DebugControl,comVecRecType,dbMAF,dbTAI93, &
       &                ComVecInfo_Init,writeComVecInfo, writeInterpsInfo,&
       &                dbStartIndex,dbEndIndex!, SnoopComVecs
  USE MLSL1Common, ONLY: Chan_R_T, Chan_R8_T, FBchans, FBnum, MBchans, MBnum, &
       WFchans, WFnum, DACSchans, DACSnum, MaxMIFs, Bandwidth, deflt_zero, R8, &
       GHzNum, BankLogical_T, BankInt_T, tau, ChanInt_T, ChanLogical_T, &
       & MaxAlts,FileNameLen, Cal_R8_T, WinMAFs, max_cal_index, &
       & space_cnts, target_cnts, limb_cnts, slimb_cnts, &
       & space_time, target_time, limb_time, slimb_time, &
       & space_weight, target_weight, slimb_weight     , &
       & space_interp, slimb_interp, target_interp,temp_interp, &
       & space_err, slimb_err, temp_err, target_err, slimb_type 


  USE L0_sci_tbls, ONLY: Sci_pkt_T
  USE EngTbls, ONLY : Eng_MAF_T
  USE Interpolation, ONLY : QuadInterpW
  USE MLSL1Config, ONLY: MIFsGHz
  USE BandTbls, ONLY: nAlts

  USE SDPToolkit, ONLY: PGS_TD_TAItoUTC


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CalWin, CalWin_T, MAFdata_T, WeightsFlags_T, Cal_R8_T, &
       &     Chan_type_T, Chi2, Tsys, Cgain


  PUBLIC :: Calibrate, InitCalibWindow, UpdateCalVectors

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(FBchans,FBnum) = "D"       ! standard filter banks
     CHARACTER(len=1) :: MB(MBchans,MBnum) = "D"       ! mid-band filter banks
     CHARACTER(len=1) :: WF(WFchans,WFnum) = "D"       ! wide filters
     CHARACTER(len=1) :: DACS(DACSchans,DACSnum) = "D" ! DACS filters
  END TYPE Chan_type_T

  !! Weights Flags type. 

  !! <whd:comment> These get set in CalibWeightsFlags.f90 when the position of
  !! the GHz switching mirror is other than what's expected for each MIF in a
  !! MAF. The `expected' settings are controlled by the Calibration section of
  !! the L1CF file (the basename of which is set via the L1CFVERSION keyword in
  !! the ENV file passed to T.sh, and the path of which is set by the inpathl1cf
  !! macro in the .macros file passed to T.sh). see the variables
  !! '{limb,space,target,discard}MIFs'. These get packed into a variable in
  !! L1Config%Calib%GHz_seq, which is then compared to what's read from
  !! telemetry (see GHz_sw_pos in CalibWeightsFlags)
  !!
  !! However, I haven't anything in the code using any field other than the
  !! %recomp_MAF field. 
  !!
  !! </whd:comment> 

  TYPE WeightsFlags_T
     INTEGER :: MAFno = 0
     LOGICAL :: recomp_MAF = .FALSE.
     LOGICAL :: recomp_S = .FALSE.
     LOGICAL :: recomp_T = .FALSE.
  END TYPE WeightsFlags_T

  !! Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (Sci_pkt_T) :: SciPkt(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     TYPE (ChanInt_T) :: LimbAltIndx, LimbAltNo
     TYPE (ChanLogical_T) :: MinCalFlag   ! Has min slimb views per MAF for cal
     TYPE (BankLogical_T) :: BankWall
     TYPE (BankInt_T) :: BankCalInd(2) ! start & end indexes for calib
     TYPE (BankInt_T) :: WallMIF       ! MIF for start of wall
     TYPE (WeightsFlags_T) :: WeightsFlags
     REAL :: MIFprecSign(0:(MaxMIFs-1))   ! Radiance precision sign per MIF
     REAL :: RnPrecSign(0:(MaxMIFs-1), 4) ! .. per MIF, per Radiometer
     INTEGER :: start_index = 0, end_index = 0  ! start/end within cal vectors
     INTEGER :: BO_stat(MIFsGHz) = 0   ! Bright Objects status (GHz FOV)
     INTEGER :: last_MIF = 0
     INTEGER :: BandSwitch(5) = 0      ! band switch positions
     LOGICAL :: CalType                ! calibration type MAF, i.e. TRUE if this MAF 
                                       ! has Calibrations in it.
     REAL :: scGeodAngle               ! Spacecraft Geod Angle (MIF 0, radians)
     TYPE (ChanLogical_T), DIMENSION(:), POINTER :: LimbAltFlag
  END TYPE MAFdata_T

  ! Moved this to MLSL1Common.f90 so as to avoid cicular dependencies
  ! with MLSL1Debug
  ! include "Calibration.f9h"

  TYPE CalWin_T
     INTEGER :: size        ! size in MAFs
     INTEGER :: current     ! current index for new data
     INTEGER :: central     ! central index to calibrate
     TYPE (MAFdata_T) :: MAFdata(WinMAFs)
     TYPE (ChanLogical_T) :: LimbAltFlag(0:(MaxMIFs-1),WinMAFs)
  END TYPE CalWin_T

  TYPE (CalWin_T), TARGET, SAVE :: CalWin
  TYPE (MAFdata_T), POINTER :: CurMAFdata

  !! Space and Target calibration vectors:
  !! Now defined in Calibration.f9h
  !! INTEGER, PARAMETER :: max_cal_index = WinMAFs * 150 - 1

  !! Moved to MLSL1Common to avoid a circular dependency with
  !! MLSL1Debug::writeRadiancesInfo

  ! TYPE (Chan_R8_T) :: space_interp(0:MaxMIFs-1)     ! Space interpolate
  ! TYPE (Chan_R8_T) :: space_err(0:MaxMIFs-1)        ! Space error
  ! TYPE (Chan_R8_T) :: slimb_interp(0:MaxMIFs-1)     ! Space/Limb interpolate
  ! TYPE (Chan_R8_T) :: slimb_err(0:MaxMIFs-1)        ! Space/Limb error
  ! TYPE (Chan_R8_T) :: target_interp(0:MaxMIFs-1)    ! Target interpolate
  ! TYPE (Chan_R8_T) :: target_err(0:MaxMIFs-1)       ! Target error
  ! TYPE (Chan_R8_T) :: temp_interp(0:MaxMIFs-1)      ! Temporary interpolate
  ! TYPE (Chan_R8_T) :: temp_err(0:MaxMIFs-1)         ! Temporary error
  ! TYPE (ChanLogical_T) :: slimb_type

  !! Counts, times, weights:

  !! moved these to MLSL1Common. This is so that I can use them in
  !! MLSL1Debug without circular dependencies.

  ! TYPE Cal_R8_T
  !    REAL(r8) :: FB(0:max_cal_index,FBchans,FBnum)
  !    REAL(r8) :: MB(0:max_cal_index,MBchans,MBnum)
  !    REAL(r8) :: WF(0:max_cal_index,WFchans,WFnum)
  !    REAL(r8) :: DACS(0:max_cal_index,DACSchans,DACSnum)
  ! END TYPE Cal_R8_T
  ! TYPE (Cal_R8_T) :: space_cnts, target_cnts, limb_cnts, slimb_cnts
  ! TYPE (Cal_R8_T) :: space_time, target_time, limb_time, slimb_time
  ! TYPE (Cal_R8_T) :: space_weight, target_weight, slimb_weight

  TYPE Cal_Int_T
     INTEGER :: FB(0:max_cal_index,FBchans,FBnum)
     INTEGER :: MB(0:max_cal_index,MBchans,MBnum)
     INTEGER :: WF(0:max_cal_index,WFchans,WFnum)
     INTEGER :: DACS(0:max_cal_index,DACSchans,DACSnum)
  END TYPE Cal_Int_T
  TYPE (Cal_Int_T) :: space_qual, target_qual, dum_qual, slimb_qual

  !! Chi square, Tsys, Cgain:

  TYPE (Chan_R_T) :: Chi2, Tsys, Cgain

  !! Important calibration private variables:

  INTEGER, SAVE :: window_MAFs, MIFsPerMAF, last_MIF
  CHARACTER(len=1) :: CalSwSeq(0:max_cal_index) = " "
  CHARACTER(len=1) :: ComVecSwSeq(0:max_cal_index) = "D"
  CHARACTER(len=1) :: CalSwSeqLS(0:max_cal_index, MaxAlts) = " "    ! Limb/Space
  CHARACTER(len=1) :: ComVecSwSeqLS(0:max_cal_index, MaxAlts) = "D" ! > min alt
  !! Generic quantities
  INTEGER  :: cal_qual(0:max_cal_index)
  REAL(r8) :: cal_weight(0:max_cal_index)
  REAL(r8) :: cal_time(0:max_cal_index)
  REAL(r8) :: errmul(0:MaxMIFs-1)

  !<whd> A really nice thing would be to know what `comVec' signifies! The best
  ! I can come up with is 'common vector', as in something which is common to
  ! anything that needs to be calibrated, regardless of whether its in the FB,
  ! MB, WF or DACS.
  ! </whd>
  

  REAL(r8) :: comVec(0:MaxMIFS-1,0:max_cal_index)

  !! MIF based, specific quantities.
  REAL(r8) :: GHz_comVec_S(0:MaxMIFs-1,0:max_cal_index)
  REAL(r8) :: GHz_comVec_T(0:MaxMIFs-1,0:max_cal_index)
  REAL(r8) :: GHz_comVec_L(0:MaxMIFs-1,0:max_cal_index,MaxAlts)
  REAL(r8) :: GHz_errmul_S(0:MaxMIFs-1)
  REAL(r8) :: GHz_errmul_T(0:MaxMIFs-1)
  REAL(r8) :: GHz_errmul_L(0:MaxMIFs-1,MaxAlts)

CONTAINS

!=============================================================================
 SUBROUTINE InitCalibWindow
!=============================================================================


    USE MLSL1Config, ONLY: L1Config
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error

    !! Initialize Calibration window data structures

    INTEGER :: i
    CHARACTER(len=23) :: asciiUTC

    window_MAFs = L1Config%Calib%CalWindow
    MIFsPerMAF = L1Config%Calib%MIFsPerMAF
    last_MIF = MIFsPerMAF - 1
    PRINT *,'WinMAFs = ',WinMAFs
    IF (window_MAFs > winMAFs) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "CalWindow greater than MAX size - resize winMAFs and " // &
            "recompile source")
    ENDIF

    !! Initialize indexes:

    CalWin%size = window_MAFs
    CalWin%current = 0        ! Indicates nothing in the window
    CalWin%central = window_MAFs / 2 + 1
    DO i = 1, window_MAFs
        CalWin%MAFdata(i)%LimbAltFlag(0:) => CalWin%LimbAltFlag(0:,i)
    ENDDO

    !! Initialize space and target weighting vectors:

    space_weight%FB = 1.0d0
    space_weight%MB = 1.0d0
    space_weight%WF = 1.0d0
    space_weight%DACS = 1.0d0
    slimb_weight%FB = 1.0d0
    slimb_weight%MB = 1.0d0
    slimb_weight%WF = 1.0d0
    slimb_weight%DACS = 1.0d0
    target_weight%FB = 1.0d0
    target_weight%MB = 1.0d0
    target_weight%WF = 1.0d0
    target_weight%DACS = 1.0d0

    IF (DebugControl%ComVecs) THEN
       comVecRecType="H "
       call ComVecInfo_Init(comVecRecType, maxMIFs,maxAlts,max_cal_index)
    ENDIF
  END SUBROUTINE InitCalibWindow

!=============================================================================
  SUBROUTINE SetComVecs
!=============================================================================
    !<whd> The only field of any _weights type that's accessed is the %FB
    !field. Why no %MB, %WF or %DACS?
    !</whd>


    USE MLSMessageModule, ONLY: MLSMessage,MLSMSG_Info

    INTEGER :: i, last_cal_index, n
    REAL(r8) :: MIFno
    CHARACTER(len=27) :: asciiUTC
    CHARACTER (LEN=FileNameLen) :: msg

    !<whd> This is for reportage, but I don't see it used anywhere else. </whd>

    IF (ANY (CalWin%MAFdata%WeightsFlags%recomp_MAF)) THEN
       PRINT *, 'win flags: ', CalWin%MAFdata%WeightsFlags
       PRINT *, 'cal type: ', CalWin%MAFdata%CalType

       write(msg,*) 'win flags: ', CalWin%MAFdata%WeightsFlags
       call MLSMessage(MLSMSG_Info,moduleName,trim(msg))
       write(msg,*) 'cal type: ', CalWin%MAFdata%CalType
       call MLSMessage(MLSMSG_Info,moduleName,trim(msg))
    ENDIF
    last_cal_index = CalWin%MAFdata(CalWin%current)%end_index

    MIFno = CalWin%MAFdata(CalWin%central)%start_index

    !<whd:debugging> MLSL1Debug
    ! pick up the MAF time of the start of this cal window
    ! and the index number of start/end
    dbStartIndex = CalWin%MAFdata(1)%start_index        
    dbEndIndex = CalWin%MAFdata(CalWin%current)%end_index
    dbMAF = CalWin%MAFdata(1)%SciPkt(0)%MAFno    
    dbTAI93 = CalWin%MAFdata(1)%SciPkt(0)%secTAI 

    n = PGS_TD_TAItoUTC (dbTAI93, asciiUTC)

    !</whd:debugging>


    !! Limb/Space: <whd> Space in Limb view? </whd>
    DO n = 1, nAlts
       IF (.NOT. (ALL (CalSwSeqLS(:,n) == ComVecSwSeqLS(:,n)))) THEN
          
          WRITE(msg,&
               & '("Doing L/S comvecs: UTC: ",a27,", MAF/index=",i4,"/",i3,", altNum: ",i2)') &
               &       asciiUTC,dbMAF,dbStartIndex,n

          PRINT *,TRIM(msg)
          CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))

          comVecRecType="LS"

          ComVecSwSeqLS(:,n) = CalSwSeqLS(:,n)   ! Save for next call

          cal_weight = slimb_weight%FB(:,1,1) ! <whd> no dependency on 'n'!</whd>

          CALL CalcComVecs (cal_qual, &            ! out
               &            cal_time, &            ! out
               &            cal_weight, &          ! in 
               &            ComVecSwSeqLS(:,n), &  ! in 
               &            "L", &                 ! in 
               &            MaxMIFs, &             ! in 
               &            last_cal_index, &      ! in 
               &            last_MIF, &            ! in 
               &            MIFno, &               ! in 
               &            comVec, &              ! out
               &            errmul)                ! out


          DO i = 0, last_MIF
             GHz_comVec_L(i,:,n) = comVec(i,:)
             GHz_errmul_L(i,n) = errmul(i)
          ENDDO

          CALL writeComVecInfo(n,dbMAF,&
               & dbtai93,CalSwSeqLS(:,n),cal_qual,cal_weight,&
               & cal_time,errmul,comVec)

       ENDIF
    ENDDO ! loop over altitudes

    IF (ALL (CalSwSeq == ComVecSwSeq)) RETURN  ! Check with current sequence

    PRINT *,'Found differences between CalSwSeq and ComVecSwSeq!'
    DO n=0,max_cal_index-1 
      IF ( CalSwSeq(n) /= ComVecSwSeq(n) ) THEN 
        PRINT *,n,CalSwSeq(n),', ',ComVecSwSeq(n)
      ENDIF
    ENDDO
    WRITE(msg,'("Doing comvecs: UTC: ",a27,", MAF/index=",i4,"/",i3)') &
         &       asciiUTC,dbMAF,dbStartIndex
    PRINT *,TRIM(msg)
    CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
    ComVecSwSeq = CalSwSeq   ! Save for next call

    !! GHz vectors:

    !! Space: ==================================================
    cal_weight = space_weight%FB(:,1,1)
    CALL CalcComVecs (cal_qual, cal_time, cal_weight, CalSwSeq, "S", &
       MaxMIFs, last_cal_index, last_MIF, MIFno, comVec, errmul)

    DO i = 0, last_MIF
       GHz_comVec_S(i,:) = comVec(i,:)
       GHz_errmul_S(i) = errmul(i)
    ENDDO

    comVecRecType="S "
    n=-1
    ! This routine does nothing unless MLSL1Debug is set!
    ! See MLSL1Debug module
    CALL writeComVecInfo(n,dbMAF,&
         &               dbtai93,CalSwSeq,cal_qual,cal_weight,&
         &               cal_time,errmul,comVec)


    !! Target ==================================================
    cal_weight = target_weight%FB(:,1,1)
    CALL CalcComVecs (cal_qual, cal_time, cal_weight, ComVecSwSeq, "T", &
       MaxMIFs, last_cal_index, last_MIF, MIFno, comVec, errmul)


    DO i = 0, last_MIF
       GHz_comVec_T(i,:) = comVec(i,:)
       GHz_errmul_T(i) = errmul(i)
    ENDDO
    comVecRecType="T "
    n=-1
    ! This routine does nothing unless MLSL1Debug is set!
    ! See MLSL1Debug module
    CALL writeComVecInfo(n,dbMAF,&
         &               dbtai93,CalSwSeq,cal_qual,cal_weight,&
         &               cal_time,errmul,comVec)

  END SUBROUTINE SetComVecs

!=============================================================================
  SUBROUTINE CalcComVecs (cal_qual, &       ! out
       &                  cal_time, &       ! out
       &                  cal_weight, &     ! in
       &                  Seq, &            ! in
       &                  cal_type, &       ! in (e.g. 'L','S','T', one 'L' could mean 'LS')
       &                  MaxMIFs, &        ! in
       &                  last_cal_index, & ! in
       &                  last_MIF, &       ! in 
       &                  MIFtime, &        ! in
       &                  comVec, &         ! out
       &                  errmul)           ! out
!=============================================================================

    INTEGER, INTENT (IN)          :: last_cal_index, last_MIF, MaxMIFs
    INTEGER, INTENT (OUT)         :: cal_qual(0:)
    REAL(r8), INTENT (OUT)        :: cal_time(0:)
    REAL(r8), INTENT (IN)         :: cal_weight(0:)
    REAL(r8), INTENT (OUT)        :: comVec(0:(MaxMIFs-1),0:last_cal_index)
    REAL(r8), INTENT (OUT)        :: errmul(0:)
    CHARACTER(len=1), INTENT (IN) :: Seq(0:), cal_type
    REAL(r8), INTENT (IN)         :: MIFtime

    INTEGER :: i, nVec, status
    REAL(r8) :: MIFno

    cal_time = -1
    cal_qual = 0
    DO i = 0, last_cal_index
       IF (Seq(i) == cal_type) THEN
          cal_time(i) = i
          cal_qual(i) = 1
       ENDIF
    ENDDO

    nVec = last_cal_index + 1

    MIFno = MIFtime

    DO i = 0, last_MIF

       CALL QuadInterpW (cal_time, cal_weight, cal_qual, MIFno, nVec, &
            comVec(i,:), errmul(i), status)

       !! Next MIF in central MAF:

       MIFno = MIFno + 1

    ENDDO

  END SUBROUTINE CalcComVecs

!!!=============================================================================
  SUBROUTINE SetCalVectors (cal_type, &   ! in
       &                    last_MIF, &   ! in
       &                    MIF_offset, & ! in
       &                    cal_cnts, &   ! out
       &                    cal_qual, &   ! out
       &                    cal_time, &   ! out
       &                    cal_altmask)  ! in (optional) Used in L/S type cals
!!!=============================================================================
    
    ! <whd>
    !
    ! Set the calibration vectors for a particular data type and
    ! MAF. Pulls data out of CurMAFdata%SciPkt(MIF)%X (X=FB,MB,WF or
    ! DACS), depending on `cal_type' and puts it into the proper
    ! location of cal_cnts.  This `unwraps' the data (which is
    ! nMAFsPerMIF x WinMAFs) and puts it into an vector that is
    ! maxMIFs*WinMIFs (maxMIFs=150, normally). The data is
    ! concatenated, so if there are (as is normally the case) 148 MIFs
    ! in all the MAFs in the calibration window, than then usable data
    ! in cal_cnts is between indices 0 and 1479.
    !
    ! </whd>

    CHARACTER(LEN=1), INTENT (IN) :: cal_type 
    INTEGER, INTENT (IN) :: last_MIF
    INTEGER, INTENT (IN) :: MIF_offset
    TYPE (Cal_R8_T), INTENT (OUT) :: cal_cnts, cal_time
    TYPE (Cal_Int_T), INTENT (OUT) :: cal_qual
    TYPE (ChanLogical_T), INTENT(IN), OPTIONAL :: cal_altmask(0:)

    INTEGER :: MIF
    TYPE (ChanLogical_T) :: altmask(0:last_MIF)

!! Initialize to not available:

    cal_time%FB(MIF_offset:MIF_offset+last_MIF,:,:) = -1  ! not available
    cal_cnts%FB(MIF_offset:MIF_offset+last_MIF,:,:) = 0
    cal_qual%FB(MIF_offset:MIF_offset+last_MIF,:,:) = 0   ! don't use for interp
    cal_time%MB(MIF_offset:MIF_offset+last_MIF,:,:) = -1  ! not available
    cal_cnts%MB(MIF_offset:MIF_offset+last_MIF,:,:) = 0
    cal_qual%MB(MIF_offset:MIF_offset+last_MIF,:,:) = 0   ! don't use for interp
    cal_time%WF(MIF_offset:MIF_offset+last_MIF,:,:) = -1  ! not available
    cal_cnts%WF(MIF_offset:MIF_offset+last_MIF,:,:) = 0
    cal_qual%WF(MIF_offset:MIF_offset+last_MIF,:,:) = 0   ! don't use for interp
    cal_time%DACS(MIF_offset:MIF_offset+last_MIF,:,:) = -1  ! not available
    cal_cnts%DACS(MIF_offset:MIF_offset+last_MIF,:,:) = 0
    cal_qual%DACS(MIF_offset:MIF_offset+last_MIF,:,:) = 0   ! don't use

!! Determine what altitude to check against, if any

    IF (PRESENT (cal_altmask)) THEN
       altmask = cal_altmask
    ELSE
       DO MIF = 0, last_MIF
          altmask(MIF)%FB = .TRUE.    ! everything
          altmask(MIF)%MB = .TRUE.    ! should
          altmask(MIF)%WF = .TRUE.    ! go
          altmask(MIF)%DACS = .TRUE.  ! through
       ENDDO
    ENDIF

    ! <whd> where the MIF has a cal, the altitutdes aren't out of range and the
    ! MIF numbers agree. </whd>

    DO MIF = 0, last_MIF

       WHERE (CurMAFdata%ChanType(MIF)%FB == cal_type .AND. altmask(MIF)%FB)
               cal_time%FB(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%FB(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%FB
          cal_qual%FB(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%MB == cal_type .AND. altmask(MIF)%MB)
          cal_time%MB(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%MB(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%MB
          cal_qual%MB(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%WF == cal_type .AND. altmask(MIF)%WF)
          cal_time%WF(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%WF(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%WF
          cal_qual%WF(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%DACS == cal_type .AND. altmask(MIF)%DACS)
          cal_time%DACS(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%DACS(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%DACS
          cal_qual%DACS(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

    ENDDO

  END SUBROUTINE SetCalVectors

!=============================================================================
  SUBROUTINE UpdateCalVectors
!=============================================================================

    !! Update the Calibration Vectors used by the interpolator. Called
    !! from SortQualify::SortAndQualify
    
    INTEGER :: windx, MIF_offset, last_MIF, start_index, end_index
    INTEGER :: bankno, channo, AltNo

    DO windx = 1, CalWin%current

       CurMAFdata => CalWin%MAFdata(windx)             ! point to current data
       last_MIF = CalWin%MAFdata(windx)%last_MIF       ! 0-based 'end-of-MAF' MIF index
       MIF_offset = CalWin%MAFdata(windx)%start_index  ! MIF offset from
                                                       ! beginning of
                                                       ! Calibration Window
       !! Update Space vectors:

       CALL SetCalVectors ("S", last_MIF, MIF_offset, space_cnts, space_qual, &
            space_time)

       !! Update Target vectors:

       CALL SetCalVectors ("T", last_MIF, MIF_offset, target_cnts, &
            target_qual, target_time)

       !! Update Limb vectors:

       CALL SetCalVectors ("L", last_MIF, MIF_offset, limb_cnts, dum_qual, &
            limb_time)

       !! Update Limb/Space vectors:

       CALL SetCalVectors ("L", last_MIF, MIF_offset, slimb_cnts, slimb_qual, &
            slimb_time, CurMAFdata%LimbAltFlag(0:last_MIF))

       !! Update Switching Sequence:

       start_index = CalWin%MAFdata(windx)%start_index
       end_index = CalWin%MAFdata(windx)%end_index

       !! <whd> These are for the non-LS type calibrations</whd>
       CalSwSeq(start_index:end_index) = &
            CalWin%MAFdata(windx)%SciPkt(0:last_MIF)%GHz_Sw_pos

       
       !! <whd> These are for the LS type calibrations, hence the need to loop
       !! over the altitudes. </whd>
       CalSwSeqLS(start_index:end_index,:) = "D"
       outer: DO AltNo = 1, nAlts
          DO bankno = 1, GHzNum
             DO channo = 1, FBchans
                IF (CurMAFdata%LimbAltIndx%FB(channo,bankno) == AltNo) THEN
                   WHERE (CurMAFdata%LimbAltFlag(0:last_MIF)%FB(channo,bankno) &
                        .AND. CurMAFdata%ChanType(0:last_MIF)%FB(channo,bankno)&
                        == "L")
                      CalSwSeqLS(start_index:end_index,AltNo) = "L"
                   ENDWHERE
                   CYCLE outer
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, MBnum
             DO channo = 1, MBchans
                IF (CurMAFdata%LimbAltIndx%MB(channo,bankno) == AltNo) THEN
                   WHERE (CurMAFdata%LimbAltFlag(0:last_MIF)%MB(channo,bankno) &
                        .AND. CurMAFdata%ChanType(0:last_MIF)%MB(channo,bankno)&
                        == "L")
                      CalSwSeqLS(start_index:end_index,AltNo) = "L"
                   ENDWHERE
                   CYCLE outer
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, WFnum
             DO channo = 1, WFchans
                IF (CurMAFdata%LimbAltIndx%WF(channo,bankno) == AltNo) THEN
                   WHERE (CurMAFdata%LimbAltFlag(0:last_MIF)%WF(channo,bankno) &
                        .AND. CurMAFdata%ChanType(0:last_MIF)%WF(channo,bankno)&
                        == "L")
                      CalSwSeqLS(start_index:end_index,AltNo) = "L"
                   ENDWHERE
                   CYCLE outer
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, DACSnum
             DO channo = 1, DACSchans
                IF (CurMAFdata%LimbAltIndx%DACS(channo,bankno) == AltNo) THEN
                   WHERE ( &
                        CurMAFdata%LimbAltFlag(0:last_MIF)%DACS(channo,bankno) &
                        .AND. &
                        CurMAFdata%ChanType(0:last_MIF)%DACS(channo,bankno) &
                        == "L")
                      CalSwSeqLS(start_index:end_index,AltNo) = "L"
                   ENDWHERE
                   CYCLE outer
                ENDIF
             ENDDO
          ENDDO
       ENDDO outer

    ENDDO

  END SUBROUTINE UpdateCalVectors

!=============================================================================
  SUBROUTINE InterpCals (calMAFs, &    ! in
       &                 nVec, &       ! in 
       &                 time, &       ! in 
       &                 cal_cnts, &   ! in 
       &                 cal_interp, & ! out
       &                 cal_err, &    ! out
       &                 GHz_comVec, & ! in 
       &                 GHz_errmul, & ! in 
       &                 BankCalInd, & ! in 
       &                 cal_index, &  ! in 
       &                 cal_time, &   ! in 
       &                 cal_weight, & ! in 
       &                 cal_qual)     ! in
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    INTEGER, INTENT(IN) :: calMAFs, nVec, cal_index(2)
    REAL(r8), INTENT(IN) :: time, GHz_errmul
    TYPE (Cal_R8_T), INTENT(IN) :: cal_cnts, cal_time, cal_weight
    TYPE (Cal_Int_T), INTENT(IN) :: cal_qual
    TYPE (Chan_R8_T),INTENT(OUT) :: cal_interp, cal_err
    REAL(r8), INTENT(IN) :: GHz_comVec(0:nVec-1)
    TYPE (BankInt_T), INTENT(IN) :: BankCalInd(2) ! start & end indexes for cals

    INTEGER :: i, j, Istat, cal1, cal2, calen
    REAL(r8) :: errmul
    REAL(r8), TARGET :: comVec(0:Max_cal_index)

    !! Previous arguments
    
    REAL(r8), SAVE, TARGET :: comVecP(0:Max_cal_index)
    REAL(r8), SAVE :: tVecP(0:Max_cal_index)
    REAL(r8), SAVE :: timeP, errmulP
    INTEGER, SAVE :: qualVecP(0:Max_cal_index)
    INTEGER, SAVE :: statusP

    REAL(r8), DIMENSION (:), POINTER :: comVecPtr

    LOGICAL :: is_same
    INTEGER, PARAMETER :: MinCalMAFs = 3  ! minimum calibration MAFs needed

! vpp orig: Set these until futher notice since comvec is calculated on a MAF basis
! before calling this routine

! <whd>
!
! To the casual observer, the logical `is_same' might suggest itself as a means
! of testing to see whether some calculation should be redone because something
! had changed. In versions of this module up to 2.3, that was the case, there
! were tests inside this routine that set `is_same' on the basis of whether
! `time' and `quality' vectors had changed. However, as of 2.3, it was
! peremptorily set to .TRUE. and has remained so since then. The only thing that
! can force calls to InterpCals to be executed is whether the length of the
! vector quantities (nVec) has changed since the last call.
!
! Revision 2.3 was committed in 2002 and, since this part of it hasn't changed
! since then, I'm assuming that it's correct, even though the effect below is to
! bypass most (all, in the case I'm working with) of the calls to InterpCals
! below.
!
! So, if you're looking for a reason why something isn't working, I'm guessing
! that looking at `is_same' is looking in the wrong direction.
!
! </whd>

    comVecP(0:nVec-1) = GHz_comVec   !! TEST!!!
    errmulP = GHz_errmul   !! TEST!!!
    is_same = .TRUE.       !! TEST!!!

    cal_interp%FB = 0.0
    cal_err%FB = 0.0
    cal_interp%MB = 0.0
    cal_err%MB = 0.0
    cal_interp%WF = 0.0
    cal_err%WF = 0.0
    cal_interp%DACS = 0.0
    cal_err%DACS = 0.0

    IF (calMAFs < MinCalMAFs) RETURN  ! Nothing can be interpolated

! Interpolate calibration values

    DO j = 1, GHzNum

       cal1 = BankCalInd(1)%FB(j)
       cal2 = BankCalInd(2)%FB(j)
       calen = cal2 - cal1 + 1

       !<whd>
       !
       ! cal[12] are the start and stop of the useable data for each of
       ! {FB,MB,WF,DACS}. Normally this will equal the start/end MIF of the
       ! current MAF unless there are BankWalls, or it could equal (0,0)
       ! However, I'm dubious about the BankWall code SortQualify::QualifyWindow
       ! (~line 1057)

       ! cal_index stores the indices of the starting MIF of the MAF before the
       ! central window and the MIF at the end of the central window MAF. This
       ! test is saying "If the start of calibration window is less than 1 MAF
       ! of the beginning of the central window, *or* it's actually in the
       ! central window, skip it. The ATB mentions ignoring MIFs close to the
       ! measurement you're calibration.
       !
       ! This test also handles the case where a bank/band is off (e.g. band
       ! 13). In that instance, cal1==cal2==0 and, therefore, cal2 <
       ! cal_index(2), so cal_interp%{FB,MB,WF,DACS} is set to 0.0
       !
       !</whd>

       IF (cal1 > cal_index(1) .OR. cal2 < cal_index(2)) THEN
          cal_interp%FB(:,j) = 0.0
          cal_err%FB(:,j) = 0.0
       ELSE
          DO i = 1, FBchans

             IF (is_same .AND. calen == nVec) THEN
                ! <whd> Nothing's changed, so just point at the last calculated comVec</whd>
                comVecPtr => comVecP
                errmul = errmulP
                Istat = statusP
             ELSE

                CALL QuadInterpW (cal_time%FB(cal1:cal2,i,j), &
                     cal_weight%FB(cal1:cal2,i,j), &
                     cal_qual%FB(cal1:cal2,i,j), time, calen, &
                     comVec(0:calen-1), errmul, Istat)

                comVecPtr => comVec

                !! Save for next call:

                tVecP = cal_time%FB(:,i,j)
                qualVecP = cal_qual%FB(:,i,j)
                timeP = time
                statusP = Istat
             ENDIF
             ! Eq D.4/D.5 of ATB
             cal_interp%FB(i,j) = SUM (comVecPtr(0:calen-1) * &
                  cal_cnts%FB(cal1:cal2,i,j))
             cal_err%FB(i,j) = errmul
          ENDDO ! loop over FB chan
       ENDIF
    ENDDO  ! GHz num

    DO j = 1, MBnum
       cal1 = BankCalInd(1)%MB(j)
       cal2 = BankCalInd(2)%MB(j)
       calen = cal2 - cal1 + 1

       IF (cal1 > cal_index(1) .OR. cal2 < cal_index(2)) THEN
          cal_interp%MB(:,j) = 0.0
          cal_err%MB(:,j) = 0.0
       ELSE
          DO i = 1, MBchans

             IF (is_same .AND. calen == nVec) THEN
                comVecPtr => comVecP
                errmul = errmulP
                Istat = statusP
             ELSE
                CALL QuadInterpW (cal_time%MB(cal1:cal2,i,j), &
                     cal_weight%MB(cal1:cal2,i,j), cal_qual%MB(cal1:cal2,i,j), &
                     time, calen, comVec(0:calen-1), errmul, Istat)

                comVecPtr => comVec

                !! Save for next call:

                tVecP = cal_time%MB(:,i,j)
                qualVecP = cal_qual%MB(:,i,j)
                timeP = time
                statusP = Istat
             ENDIF
             ! Eq D.5 of ATB (D.4 is implemented in QuadInterpW)
             cal_interp%MB(i,j) = SUM (comVecPtr(0:calen-1) * &
                  cal_cnts%MB(cal1:cal2,i,j))
             cal_err%MB(i,j) = errmul
          ENDDO ! MB chan
       ENDIF
    ENDDO  ! MB num

    DO j = 1, WFnum
       cal1 = BankCalInd(1)%WF(j)
       cal2 = BankCalInd(2)%WF(j)
       calen = cal2 - cal1 + 1

       IF (cal1 > cal_index(1) .OR. cal2 < cal_index(2)) THEN
          cal_interp%WF(:,j) = 0.0
          cal_err%WF(:,j) = 0.0
       ELSE
          DO i = 1, WFchans

             IF (is_same .AND. calen == nVec) THEN
                comVecPtr => comVecP
                errmul = errmulP
                Istat = statusP
             ELSE
                CALL QuadInterpW (cal_time%WF(cal1:cal2,i,j), &
                     cal_weight%WF(cal1:cal2,i,j), cal_qual%WF(cal1:cal2,i,j), &
                     time, calen, comVec(0:calen-1), errmul, Istat)

                comVecPtr => comVec

                !! Save for next call:

                tVecP = cal_time%WF(:,i,j)
                qualVecP = cal_qual%WF(:,i,j)
                timeP = time
                statusP = Istat
             ENDIF
             cal_interp%WF(i,j) = SUM (comVecPtr(0:calen-1) * &
                  cal_cnts%WF(cal1:cal2,i,j))
             cal_err%WF(i,j) = errmul
          ENDDO ! Loop over WF channels
       ENDIF
    ENDDO  ! WF num

    IF (L1Config%Calib%CalibDACS) THEN

       DO j = 1, DACSnum
          cal1 = BankCalInd(1)%DACS(j)
          cal2 = BankCalInd(2)%DACS(j)
          calen = cal2 - cal1 + 1
 
          IF (cal1 > cal_index(1) .OR. cal2 < cal_index(2)) THEN
             cal_interp%DACS(:,j) = 0.0
             cal_err%DACS(:,j) = 0.0
          ELSE
             DO i = 1, DACSchans

                IF (is_same .AND. calen == nVec) THEN
                   comVecPtr => comVecP
                   errmul = errmulP
                   Istat = statusP
                ELSE
                   CALL QuadInterpW (cal_time%DACS(cal1:cal2,i,j), &
                        cal_weight%DACS(cal1:cal2,i,j), &
                        cal_qual%DACS(cal1:cal2,i,j), time, calen, &
                        comVec(0:calen-1), errmul, Istat)

                   comVecPtr => comVec

                   !! Save for next call:

                   tVecP = cal_time%DACS(:,i,j)
                   qualVecP = cal_qual%DACS(:,i,j)
                   timeP = time
                   statusP = Istat
                ENDIF
                ! Eq D.4/D.5 of ATB
                cal_interp%DACS(i,j) = SUM (comVecPtr(0:calen-1) * &
                     cal_cnts%DACS(cal1:cal2,i,j))
                cal_err%DACS(i,j) = errmul
             ENDDO ! Loop over channel
          ENDIF
       ENDDO ! Loop over DACS number

    ENDIF

  END SUBROUTINE InterpCals

!=============================================================================
  SUBROUTINE ChiSquare (start_index, end_index, space_counts, space_interp, &
       nlast)
!=============================================================================

    ! <whd> 
    ! Calculates the quantity that's eventually written to the L1BDIAG
    ! file as /Chi2 {FB,MB,WF,DACS}'
    ! </whd>

    INTEGER :: start_index, end_index, nlast
    TYPE (Cal_R8_T) :: space_counts
    TYPE (Chan_R8_T) :: space_interp(0:nlast)
    INTEGER :: i, j, nspace
    INTEGER :: nvec(0:nlast)
    REAL(r8) :: difspace(0:nlast), difzero(0:nlast)
    REAL(r8) :: SumDifS2, SumDifZ2
    INTEGER, PARAMETER :: minmafs = 6

    chi2%FB = 0.0      ! initial value
    DO j = 1, FBnum
       DO i = 1, FBchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%FB(start_index:end_index,i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%FB(start_index:end_index,i,j) - &
                  space_interp%FB(i,j)
             difzero = space_counts%FB(start_index:end_index,i,j) - &
                  deflt_zero%FB(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%FB(i,j) = &
               & (SumDifS2 / nspace - (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
                     / (bandwidth%FB(i,j) * tau))
                chi2%FB(i,j) = chi2%FB(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    chi2%MB = 0.0      ! initial value
    DO j = 1, MBnum
       DO i = 1, MBchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%MB(start_index:end_index,i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%MB(start_index:end_index,i,j) - &
                  space_interp%MB(i,j)
             difzero = space_counts%MB(start_index:end_index,i,j) - &
                  deflt_zero%MB(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%MB(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) / ((SumDifZ2 / nspace) &
                     / (bandwidth%MB(i,j) * tau))
                chi2%MB(i,j) = chi2%MB(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    chi2%WF = 0.0      ! initial value
    DO j = 1, WFnum
       DO i = 1, WFchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%WF(start_index:end_index,i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%WF(start_index:end_index,i,j) - &
                  space_interp%WF(i,j)
             difzero = space_counts%WF(start_index:end_index,i,j) - &
                  deflt_zero%WF(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%WF(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
                     / (bandwidth%WF(i,j) * tau))
                chi2%WF(i,j) = chi2%WF(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    chi2%DACS = 0.0      ! initial value
    DO j = 1, DACSnum
       DO i = 1, DACSchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%DACS(start_index:end_index,i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%DACS(start_index:end_index,i,j) - &
                  space_interp%DACS(i,j)
             difzero = space_counts%DACS(start_index:end_index,i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%DACS(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
                     / (bandwidth%DACS(i,j) * tau))
                chi2%DACS(i,j) = chi2%DACS(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ChiSquare

!=============================================================================
  SUBROUTINE Calibrate
!=============================================================================

    USE MLSL1Rad, ONLY: UpdateRadSignals
    USE MLSL1Config, ONLY: L1Config

!! Calibrate the science data

    INTEGER :: time_index, start_index, end_index, windex
    INTEGER :: AltNo, bankno, channo, nVec, cal_index(2), MIF_index, calMAFs
    REAL(r8) :: time, secs, tai93
    REAL(r8), SAVE :: oldsecs 
    data oldsecs/0.0/

    CHARACTER(len=8) :: date
    CHARACTER (len=10) :: dtime
    CHARACTER (len=5) :: zone
    INTEGER :: values(8)
    TYPE (ChanInt_T), POINTER :: LimbAltIndx          ! Channel Limb indexes
    LOGICAL :: do_slimb = .FALSE.

    PRINT *, 'calibrating...'

    ! print out some execution time information. 
    CALL DATE_AND_TIME (date, dtime, zone, values)
    secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
    IF (oldsecs == 0.0) THEN 
       print *,'First time through! No timing information available'
    ELSE
       PRINT *, "Time between calls to this routine: ", (secs-oldsecs)
    ENDIF
    oldsecs = secs

    nVec = CalWin%MAFdata(CalWin%size)%end_index + 1
    windex = CalWin%central ! 'central' MAF in WinMAFs worth of data.

    ! <whd:comment> 
    !
    ! {start,end}_index marks the beginning/ending of the `central' MAF in the
    ! WinMAFs list of MAFs in the calibration window.
    !
    ! cal_index stores the indices of the beginning of the MAF before
    ! start_index and to the end of the central MAF
    !
    ! </whd:comment> 

    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    tai93=CalWin%MAFdata(windex)%SciPkt(0)%secTAI     ! start of cal window TAI93 time 
    cal_index(1) = start_index - CalWin%MAFdata(windex-1)%EMAF%MIFsPerMAF
    cal_index(2) = start_index + CalWin%MAFdata(windex)%EMAF%MIFsPerMAF

    !PRINT *,'{start,end}_index = ',start_index,end_index
    !PRINT *,'cal_index = ',cal_index
    ! Indices of MIFs with good altitudes.(???)
    LimbAltIndx => CalWin%MAFdata(windex)%LimbAltIndx
    ! <whd:comment>
    !
    ! The 'signal' refered to here is a complicated user type:
    ! MLSSignalNomenclature::MLSSignal_T which stores information about what
    ! data is moving through what bit of electronics based on the value of
    ! BandSwitch. The 'signal' doesn't have any telemetry data in it.
    !
    ! </whd:comment>
    CALL UpdateRadSignals (CalWin%MAFdata(windex)%BandSwitch)
    ! 
    CALL SetComVecs

    calMAFs = COUNT (CalWin%MAFdata%CalType)

    DO time_index = start_index, end_index  ! for every MIF in the central MAF

       time = time_index                     ! "Time" from start of cal window
       MIF_index = time_index - start_index  ! MIF # within the central MAF

       ! Space cals:

       ! <whd> Everything in this call *except* space_{interp,err} are
       ! inputs. Same, Same for next two calls.
       ! </whd>

       CALL InterpCals (calMAFs, nVec, time, space_cnts, &
            space_interp(MIF_index), space_err(MIF_index), &
            GHz_comVec_S(MIF_index,:), GHz_errmul_S(MIF_index), &
            CalWin%MAFdata(windex)%BankCalInd, cal_index, space_time, &
            space_weight, space_qual)

       ! Target cals:

       CALL InterpCals (calMAFs, nVec, time, target_cnts, &
            target_interp(MIF_index), target_err(MIF_index), &
            GHz_comVec_T(MIF_index,:), GHz_errmul_T(MIF_index), &
            CalWin%MAFdata(windex)%BankCalInd, cal_index, target_time, &
            target_weight, target_qual)

       ! Limb/Space cals:

       DO AltNo = 1, nAlts

          CALL InterpCals (calMAFs, nVec, time, slimb_cnts, &
               temp_interp(MIF_index), temp_err(MIF_index), &
               GHz_comVec_L(MIF_index,:,AltNo), GHz_errmul_L(MIF_index,AltNo), &
               CalWin%MAFdata(windex)%BankCalInd, cal_index, slimb_time, &
               slimb_weight, slimb_qual)

          ! Save appropriate channel interpolates:

          DO bankno = 1, GHzNum
             DO channo = 1, FBchans
                IF (LimbAltIndx%FB(channo,bankno) == AltNo) THEN
                   slimb_interp(MIF_index)%FB(channo,bankno) = &
                        temp_interp(MIF_index)%FB(channo,bankno)
                   slimb_err(MIF_index)%FB(channo,bankno) = &
                        temp_err(MIF_index)%FB(channo,bankno)
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, MBnum
             DO channo = 1, MBchans
                IF (LimbAltIndx%MB(channo,bankno) == AltNo) THEN
                   slimb_interp(MIF_index)%MB(channo,bankno) = &
                        temp_interp(MIF_index)%MB(channo,bankno)
                   slimb_err(MIF_index)%MB(channo,bankno) = &
                        temp_err(MIF_index)%MB(channo,bankno)
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, WFnum
             DO channo = 1, WFchans
                IF (LimbAltIndx%WF(channo,bankno) == AltNo) THEN
                   slimb_interp(MIF_index)%WF(channo,bankno) = &
                        temp_interp(MIF_index)%WF(channo,bankno)
                   slimb_err(MIF_index)%WF(channo,bankno) = &
                        temp_err(MIF_index)%WF(channo,bankno)
                ENDIF
             ENDDO
          ENDDO
          DO bankno = 1, DACSnum
             DO channo = 1, DACSchans
                IF (LimbAltIndx%DACS(channo,bankno) == AltNo) THEN
                   slimb_interp(MIF_index)%DACS(channo,bankno) = &
                        temp_interp(MIF_index)%DACS(channo,bankno)
                   slimb_err(MIF_index)%DACS(channo,bankno) = &
                        temp_err(MIF_index)%DACS(channo,bankno)
                ENDIF
             ENDDO
          ENDDO
       ENDDO ! Loop over unique altitudes

    ENDDO ! Loop over MIFs in central MAF

    IF (debugControl%Interps) THEN 
      CALL writeInterpsInfo(maxMIFs, &
           &  tai93, &
           &  space_interp%FB(1,1),&
           &  space_err%FB(1,1),&
           &  slimb_interp%FB(1,1),&
           &  slimb_err%FB(1,1),&
           &  target_interp%FB(1,1),&
           &  target_err%FB(1,1))
    ENDIF

    slimb_type%FB = .FALSE.
    slimb_type%MB = .FALSE.
    slimb_type%WF = .FALSE.
    slimb_type%DACS = .FALSE.


! Mark good "slimb" channels (if requested):

    IF (L1Config%Calib%do_slimb) THEN   ! slimb_type
       WHERE (LimbAltIndx%FB /= 0 .AND. CalWin%MAFdata(windex)%MinCalFlag%FB)
          slimb_type%FB = .TRUE.
       ENDWHERE
       WHERE (LimbAltIndx%MB /= 0 .AND. CalWin%MAFdata(windex)%MinCalFlag%MB)
          slimb_type%MB = .TRUE.
       ENDWHERE
       WHERE (LimbAltIndx%WF /= 0 .AND. CalWin%MAFdata(windex)%MinCalFlag%WF)
          slimb_type%WF = .TRUE.
       ENDWHERE
       WHERE (LimbAltIndx%DACS /= 0 .AND. &
            CalWin%MAFdata(windex)%MinCalFlag%DACS)
          slimb_type%DACS = .TRUE.
       ENDWHERE

    ENDIF

    CALL ChiSquare (start_index, end_index, space_cnts, &
         space_interp(0:end_index-start_index), (end_index-start_index))

PRINT *, 'end calibrating...'

  END SUBROUTINE Calibrate

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE Calibration
!=============================================================================

! $Log$
! Revision 2.27  2024/10/10 20:17:17  pwagner
! Add new component RnPrecSign to MAFdata_T
!
! Revision 2.26  2018/04/09 22:12:58  whdaffer
! Mostly documentation, and some reportage
!
! Revision 2.25  2016/03/18 19:07:22  whdaffer
! Took out extraneous ',' that NAG complained about
!
! Revision 2.24  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.23.2.2  2016/03/14 19:51:24  whdaffer
! Most of the work is to eliminate cicular references between MLSL1Debug
! and Calibration. To resolve this, I've moved the inclusion of
! Calibration.f9h and the definition of some 10 variables from
! Calibration.f90 to MLSL1Common.f90. Radiances, MLSL1Debug and
! Calibration will get those types, variable from MLSL1Common. Also, use
! machines.f90 to get the definition of usleep used in SnoopMLSL1
!
! Revision 2.23.2.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.23  2015/04/27 20:27:56  whdaffer
! Removed obsolete preprocessor commands. Added a bit of documentation
! about Calibration.f90
!
! Revision 2.22  2015/04/27 20:25:18  whdaffer
! Moved print WinMAFs statement to InitCalibWindow
!
! Revision 2.21  2015/04/23 17:46:27  whdaffer
! removed Makefile
!
! Revision 2.20  2015/01/13 18:37:37  pwagner
! Changed lower bounds on pointer to match LimbAltFlag
!
! Revision 2.19  2007/02/09 15:02:45  perun
! Do slimb calibration only if requested.
!
! Revision 2.18  2006/09/26 16:01:05  perun
! Add DACS Chi2 calculation
!
! Revision 2.17  2006/06/14 13:44:24  perun
! Add Spacecraft Geod Angle
!
! Revision 2.16  2006/03/24 15:07:20  perun
! Add Space in Limb calibration based on limb altitude
!
! Revision 2.15  2005/12/06 19:22:30  perun
! Removed BrightObjest_T and added BO_stat to MAFdata_T
!
! Revision 2.14  2005/10/10 14:27:32  perun
! Add CalType for each MAF in CalWin and calibrate based on number of calMAFs
!
! Revision 2.13  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.12  2004/11/10 15:38:57  perun
! Moved BrightObjects_T to this routine; added MIFprecSign flag defintion
!
! Revision 2.11  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.10  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.9  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.8  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.7  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.6  2003/01/31 18:13:33  perun
! Version 1.1 commit
!
! Revision 2.4  2002/08/06 20:43:45  perun
! Set all calibration to precomputed until further notice.
!
! Revision 2.3  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.2  2001/09/10 16:16:08  perun
! Changed ALLOCATABLE component to POINTER
!
! Revision 2.1  2001/02/23 18:50:29  perun
! Version 0.5 commit
!
