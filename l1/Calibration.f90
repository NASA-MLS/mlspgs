! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Calibration ! Calibration data and routines
!=============================================================================

  USE MLSL1Common, ONLY: Chan_R_T, Chan_R8_T, FBchans, FBnum, MBchans, MBnum, &
       WFchans, WFnum, DACSchans, DACSnum, MaxMIFs, Bandwidth, deflt_zero, R8, &
       BankLogical_T, BankInt_T, tau
  USE L0_sci_tbls, ONLY: Sci_pkt_T
  USE EngTbls, ONLY : Eng_MAF_T
  USE Interpolation, ONLY : QuadInterpW

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CalWin, CalWin_T, MAFdata_T, WeightsFlags_T, Cal_R8_T, &
       Chan_type_T, BrightObjects_T
  PUBLIC :: limb_cnts, space_interp, target_interp, space_err, target_err, &
       Chi2, Tsys, Cgain

  PUBLIC :: Calibrate, InitCalibWindow, UpdateCalVectors

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(FBchans,FBnum)          ! standard filter banks
     CHARACTER(len=1) :: MB(MBchans,MBnum)          ! mid-band filter banks
     CHARACTER(len=1) :: WF(WFchans,WFnum)          ! wide filters
     CHARACTER(len=1) :: DACS(DACSchans,DACSnum)    ! DACS filters
  END TYPE Chan_type_T

  !! Weights Flags type

  TYPE WeightsFlags_T
     INTEGER :: MAFno
     LOGICAL :: recomp_MAF, recomp_S, recomp_T
  END TYPE WeightsFlags_T

  !! Bright Objects type

  TYPE BrightObjects_T
     LOGICAL :: MoonInFOV(0:MaxMIFs-1)
     LOGICAL :: VenusInFOV(0:MaxMIFs-1)
  END TYPE BrightObjects_T

  !! Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (Sci_pkt_T) :: SciPkt(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     TYPE (BankLogical_T) :: BankWall
     TYPE (BankLogical_T) :: Nominal   ! nominal switching flag
     TYPE (BankInt_T) :: BankCalInd(2) ! start & end indexes for calib
     TYPE (BankInt_T) :: WallMIF       ! MIF for start of wall
     TYPE (BrightObjects_T) :: LimbView, SpaceView ! Bright Objects in FOV flags
     TYPE (WeightsFlags_T) :: WeightsFlags
     INTEGER :: start_index, end_index  ! start & end indexes within cal vectors
     INTEGER :: last_MIF
     INTEGER :: BandSwitch(5)           ! band switch positions
  END TYPE MAFdata_T

  !! Calibration window:

  INTEGER, PARAMETER :: WinMAFs = 6     ! current window size in MAFs

  TYPE CalWin_T
     INTEGER :: size        ! size in MAFs
     INTEGER :: current     ! current index for new data
     INTEGER :: central     ! central index to calibrate
     TYPE (MAFdata_T) :: MAFdata(WinMAFs)
  END TYPE CalWin_T

  TYPE (CalWin_T), TARGET :: CalWin
  TYPE (MAFdata_T), POINTER :: CurMAFdata

  !! Space and Target calibration vectors:

  INTEGER, PARAMETER :: max_cal_index = WinMAFs * 150 - 1

  TYPE (Chan_R8_T) :: space_interp(0:MaxMIFS-1)     ! Space interpolate
  TYPE (Chan_R8_T) :: space_err(0:MaxMIFS-1)        ! Space error
  TYPE (Chan_R8_T) :: target_interp(0:MaxMIFS-1)    ! Target interpolate
  TYPE (Chan_R8_T) :: target_err(0:MaxMIFS-1)       ! Target error

  !! Counts, times, weights:

  TYPE Cal_R8_T
     REAL(r8) :: FB(0:max_cal_index,FBchans,FBnum)
     REAL(r8) :: MB(0:max_cal_index,MBchans,MBnum)
     REAL(r8) :: WF(0:max_cal_index,WFchans,WFnum)
     REAL(r8) :: DACS(0:max_cal_index,DACSchans,DACSnum)
  END TYPE Cal_R8_T
  TYPE (Cal_R8_T) :: space_cnts, target_cnts, limb_cnts
  TYPE (Cal_R8_T) :: space_time, target_time, limb_time
  TYPE (Cal_R8_T) :: space_weight, target_weight

  TYPE Cal_Int_T
     INTEGER :: FB(0:max_cal_index,FBchans,FBnum)
     INTEGER :: MB(0:max_cal_index,MBchans,MBnum)
     INTEGER :: WF(0:max_cal_index,WFchans,WFnum)
     INTEGER :: DACS(0:max_cal_index,DACSchans,DACSnum)
  END TYPE Cal_Int_T
  TYPE (Cal_Int_T) :: space_qual, target_qual, dum_qual

  !! Chi square, Tsys, Cgain:

  TYPE (Chan_R_T) :: Chi2, Tsys, Cgain

  !! Important calibration private variables:

  INTEGER, SAVE :: window_MAFs, MIFsPerMAF, last_MIF
  CHARACTER(len=1) :: CalSwSeq(0:max_cal_index) = " "
  CHARACTER(len=1) :: ComVecSwSeq(0:max_cal_index) = "D"
  INTEGER :: cal_qual(0:max_cal_index)
  REAL(r8) :: cal_weight(0:max_cal_index), cal_time(0:max_cal_index)
  REAL(r8) :: comVec(0:MaxMIFS-1,0:max_cal_index)
  REAL(r8) :: GHz_comVec_S(0:MaxMIFS-1,0:max_cal_index)
  REAL(r8) :: GHz_comVec_T(0:MaxMIFS-1,0:max_cal_index)
  REAL(r8) :: errmul(0:MaxMIFS-1)
  REAL(r8) :: GHz_errmul_S(0:MaxMIFS-1), GHz_errmul_T(0:MaxMIFS-1)

CONTAINS

!=============================================================================
 SUBROUTINE InitCalibWindow
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    !! Initialize Calibration window data structures

    INTEGER :: i

    window_MAFs = L1Config%Calib%CalWindow
    MIFsPerMAF = L1Config%Calib%MIFsPerMAF
    last_MIF = MIFsPerMAF - 1

    !! Initialize indexes:

    CalWin%size = window_MAFs
    CalWin%current = 0        ! Indicates nothing in the window
    CalWin%central = window_MAFs / 2 + 1

    !! Initialize nominal data flags

    DO i = 1, window_MAFs
       CalWin%MAFdata(i)%nominal%FB = .TRUE.
       CalWin%MAFdata(i)%nominal%MB = .TRUE.
       CalWin%MAFdata(i)%nominal%WF = .TRUE.
       CalWin%MAFdata(i)%nominal%DACS = .TRUE.
    ENDDO

    !! Initialize space and target weighting vectors:

    space_weight%FB = 1.0d0
    space_weight%MB = 1.0d0
    space_weight%WF = 1.0d0
    space_weight%DACS = 1.0d0
    target_weight%FB = 1.0d0
    target_weight%MB = 1.0d0
    target_weight%WF = 1.0d0
    target_weight%DACS = 1.0d0

  END SUBROUTINE InitCalibWindow

!=============================================================================
  SUBROUTINE SetComVecs
!=============================================================================

    INTEGER :: i, last_cal_index
    REAL(r8) :: MIFno
    LOGICAL, SAVE :: done = .FALSE.

    IF (ANY (calwin%mafdata%weightsflags%recomp_maf)) THEN
       PRINT *, 'win flags: ', CalWin%MAFdata%WeightsFlags
    ENDIF
    last_cal_index = CalWin%MAFdata(CalWin%current)%end_index

    done = ALL (CalSwSeq == ComVecSwSeq)  ! Check with current sequence

    IF (done) RETURN   ! already done

    PRINT *, 'doing comvecs...'

    ComVecSwSeq = CalSwSeq   ! Save for next call

    MIFno = CalWin%MAFdata(CalWin%central)%start_index

    !! GHz vectors:

    !! Space:

    cal_weight = space_weight%FB(:,1,1)

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, ComVecSwSeq, "S", &
       MaxMIFs, last_cal_index, last_MIF, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       GHz_comVec_S(i,:) = comVec(i,:)
       GHz_errmul_S(i) = errmul(i)
    ENDDO

    !! Target:

    cal_weight = target_weight%FB(:,1,1)

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, ComVecSwSeq, "T", &
       MaxMIFs, last_cal_index, last_MIF, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       GHz_comVec_T(i,:) = comVec(i,:)
       GHz_errmul_T(i) = errmul(i)
    ENDDO

  END SUBROUTINE SetComVecs

!=============================================================================
  SUBROUTINE CalcComVecs (cal_qual, cal_time, cal_weight, Seq, cal_type, &
       MaxMIFs, last_cal_index, last_MIF, MIFtime, comVec, errmul)
!=============================================================================

    INTEGER :: last_cal_index, last_MIF, MaxMIFs
    INTEGER :: cal_qual(0:)
    REAL(r8) :: cal_time(0:)
    REAL(r8) :: cal_weight(0:)
    REAL(r8) :: comVec(0:(MaxMIFs-1),0:last_cal_index)
    REAL(r8) :: errmul(0:)
    CHARACTER(len=1) :: Seq(0:), cal_type
    REAL(r8) :: MIFtime

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

!=============================================================================
  SUBROUTINE SetCalVectors (cal_type, last_MIF, MIF_offset, cal_cnts, &
       cal_qual, cal_time)
!=============================================================================

!! Set the calibration vectors for a particular data type

    CHARACTER(LEN=1), INTENT (IN) :: cal_type 
    INTEGER, INTENT (IN) :: last_MIF
    INTEGER, INTENT (IN) :: MIF_offset
    TYPE (Cal_R8_T), INTENT (OUT) :: cal_cnts, cal_time
    TYPE (Cal_Int_T), INTENT (OUT) :: cal_qual

    INTEGER :: MIF

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

    DO MIF = 0, last_MIF

       WHERE (CurMAFdata%ChanType(MIF)%FB == cal_type)
          cal_time%FB(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%FB(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%FB
          cal_qual%FB(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%MB == cal_type)
          cal_time%MB(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%MB(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%MB
          cal_qual%MB(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%WF == cal_type)
          cal_time%WF(MIF_offset+MIF,:,:) = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cal_cnts%WF(MIF_offset+MIF,:,:) = CurMAFdata%SciPkt(MIF)%WF
          cal_qual%WF(MIF_offset+MIF,:,:) = 1    ! use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%DACS == cal_type)
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

!! Update the Calibration Vectors used by the interpolator

    INTEGER :: windx, MIF_offset, last_MIF, start_index, end_index

    DO windx = 1, CalWin%current

       CurMAFdata => CalWin%MAFdata(windx)             ! point to current data
       last_MIF = CalWin%MAFdata(windx)%last_MIF
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

       !! Update Switching Sequence:

       start_index = CalWin%MAFdata(windx)%start_index
       end_index = CalWin%MAFdata(windx)%end_index
       CalSwSeq(start_index:end_index) = &
            CalWin%MAFdata(windx)%SciPkt(0:last_MIF)%GHz_Sw_pos

    ENDDO

  END SUBROUTINE UpdateCalVectors

!=============================================================================
  SUBROUTINE InterpCals (nVec, time, cal_cnts, cal_interp, cal_err, &
       GHz_comVec, GHz_errmul, BankCalInd, cal_index, cal_time, cal_weight, &
       cal_qual)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    INTEGER :: nVec, cal_index(2)
    REAL(r8) :: time, GHz_errmul
    TYPE (Cal_R8_T) :: cal_cnts, cal_time, cal_weight
    TYPE (Cal_Int_T) :: cal_qual
    TYPE (Chan_R8_T) :: cal_interp, cal_err
    REAL(r8) :: GHz_comVec(0:nVec-1)
    TYPE (BankInt_T) :: BankCalInd(2) ! start & end indexes for calib

    INTEGER :: i, j, Istat, cal1, cal2, calen, calMAFs
    REAL(r8) :: errmul

    !! Previous arguments
    
    REAL(r8), SAVE :: tVecP(0:Max_cal_index), comVecP(0:Max_cal_index), &
         comVec(0:Max_cal_index)
    REAL(r8), SAVE :: timeP, errmulP
    INTEGER, SAVE :: qualVecP(0:Max_cal_index)
    INTEGER, SAVE :: statusP

    LOGICAL :: is_same
    INTEGER, PARAMETER :: MinCalMAFs = 3  ! minimum calibration MAFs needed

comVecP(0:nVec-1) = GHz_comVec   !! TEST!!!
errmulP = GHz_errmul   !! TEST!!!
is_same = .true.       !! TEST!!!

! Interpolate calibration values

    DO j = 1, (FBnum-5)   ! Don't need to do THz here!!!
       cal1 = BankCalInd(1)%FB(j)
       cal2 = BankCalInd(2)%FB(j)
       calen = cal2 - cal1 + 1
       calMAFs = calen / L1Config%Calib%MIFsPerMAF

       IF (calMAFs < MinCalMAFs .OR. cal1 > cal_index(1) .OR. &
            cal2 < cal_index(2)) THEN
! if (calMAFs > MinCalMAFs) print *, 'calMAFs: ', calMAFs
          cal_interp%FB(:,j) = 0.0
          cal_err%FB(:,j) = 0.0
       ELSE
          DO i = 1, FBchans

             IF (is_same .AND. calen == nVec) THEN
                comVec = comVecP
                errmul = errmulP
                Istat = statusP
             ELSE

                CALL QuadInterpW (cal_time%FB(cal1:cal2,i,j), &
                     cal_weight%FB(cal1:cal2,i,j), &
                     cal_qual%FB(cal1:cal2,i,j), time, calen, &
                     comVec(0:calen-1), errmul, Istat)

                !! Save for next call:

                tVecP = cal_time%FB(:,i,j)
                qualVecP = cal_qual%FB(:,i,j)
                timeP = time
                statusP = Istat
             ENDIF

             cal_interp%FB(i,j) = &
                  SUM (comVec(0:calen-1) * cal_cnts%FB(cal1:cal2,i,j))
             cal_err%FB(i,j) = errmul
          ENDDO
       ENDIF
    ENDDO

    DO j = 1, MBnum
       cal1 = BankCalInd(1)%MB(j)
       cal2 = BankCalInd(2)%MB(j)
       calen = cal2 - cal1 + 1
       DO i = 1, MBchans

          IF (is_same .AND. calen == nVec) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time%MB(cal1:cal2,i,j), &
                  cal_weight%MB(cal1:cal2,i,j), cal_qual%MB(cal1:cal2,i,j), &
                  time, calen, comVec(0:calen-1), errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%MB(:,i,j)
             qualVecP = cal_qual%MB(:,i,j)
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%MB(i,j) = &
               SUM (comVec(0:calen-1) * cal_cnts%MB(cal1:cal2,i,j))
          cal_err%MB(i,j) = errmul
       ENDDO
    ENDDO

    DO j = 1, WFnum
       cal1 = BankCalInd(1)%WF(j)
       cal2 = BankCalInd(2)%WF(j)
       calen = cal2 - cal1 + 1
       DO i = 1, WFchans

          IF (is_same .AND. calen == nVec) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time%WF(cal1:cal2,i,j), &
                  cal_weight%WF(cal1:cal2,i,j), cal_qual%WF(cal1:cal2,i,j), &
                  time, calen, comVec(0:calen-1), errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%WF(:,i,j)
             qualVecP = cal_qual%WF(:,i,j)
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%WF(i,j) = &
               SUM (comVec(0:calen-1) * cal_cnts%WF(cal1:cal2,i,j))
          cal_err%WF(i,j) = errmul
       ENDDO
    ENDDO

    IF (L1Config%Calib%CalibDACS) THEN

       DO j = 1, DACSnum
          cal1 = BankCalInd(1)%DACS(j)
          cal2 = BankCalInd(2)%DACS(j)
          calen = cal2 - cal1 + 1
          DO i = 1, DACSchans

             IF (is_same .AND. calen == nVec) THEN
                comVec = comVecP
                errmul = errmulP
                Istat = statusP
             ELSE
                CALL QuadInterpW (cal_time%DACS(cal1:cal2,i,j), &
                     cal_weight%DACS(cal1:cal2,i,j), &
                     cal_qual%DACS(cal1:cal2,i,j), time, calen, &
                     comVec(0:calen-1), errmul, Istat)

                !! Save for next call:

                tVecP = cal_time%DACS(:,i,j)
                qualVecP = cal_qual%DACS(:,i,j)
                timeP = time
                statusP = Istat
             ENDIF
             cal_interp%DACS(i,j) = &
                  SUM (comVec(0:calen-1) * cal_cnts%DACS(cal1:cal2,i,j))
             cal_err%DACS(i,j) = errmul
          ENDDO
       ENDDO

    ENDIF

  END SUBROUTINE InterpCals

!=============================================================================
  SUBROUTINE ChiSquare (start_index, end_index, space_counts, space_interp, &
       nlast)
!=============================================================================

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
                chi2%FB(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
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

  END SUBROUTINE ChiSquare

!=============================================================================
  SUBROUTINE Calibrate
!=============================================================================

    USE MLSL1Rad, ONLY: UpdateRadSignals

!! Calibrate the science data

    INTEGER :: time_index, start_index, end_index, windex
    INTEGER :: nVec, cal_index(2), MIF_index
    REAL(r8) :: time, secs, oldsecs

    CHARACTER(len=8) :: date
    CHARACTER (len=10) :: dtime
    CHARACTER (len=5) :: zone
    INTEGER :: values(8)

PRINT *, 'calibrating...'

CALL DATE_AND_TIME (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
oldsecs = secs

    nVec = CalWin%MAFdata(CalWin%size)%end_index + 1
    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    cal_index(1) = start_index - CalWin%MAFdata(windex-1)%EMAF%MIFsPerMAF
    cal_index(2) = start_index + CalWin%MAFdata(windex)%EMAF%MIFsPerMAF

    CALL UpdateRadSignals (CalWin%MAFdata(windex)%BandSwitch)

    CALL SetComVecs

    DO time_index = start_index, end_index  ! for every MIF in the MAF

CALL DATE_AND_TIME (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
if (time_index == end_index) then
print *, "Time dif: ", (secs-oldsecs)
else if (time_index == start_index) then
oldsecs = secs
endif
       time = time_index                     ! "Time" from start of cal window
       MIF_index = time_index - start_index  ! MIF # within the central MAF

       ! Space cals:

       CALL InterpCals (nVec, time, space_cnts, space_interp(MIF_index), &
            space_err(MIF_index), GHz_comVec_S(MIF_index,:), &
            GHz_errmul_S(MIF_index), CalWin%MAFdata(windex)%BankCalInd, &
            cal_index, space_time, space_weight, space_qual)

       ! Target cals:

       CALL InterpCals (nVec, time, target_cnts, target_interp(MIF_index), &
            target_err(MIF_index), GHz_comVec_T(MIF_index,:), &
            GHz_errmul_T(MIF_index), CalWin%MAFdata(windex)%BankCalInd, &
            cal_index, target_time, target_weight, target_qual)

    ENDDO

    CALL ChiSquare (start_index, end_index, space_cnts, &
         space_interp(0:end_index-start_index), (end_index-start_index))

PRINT *, 'end calibrating...'

  END SUBROUTINE Calibrate

!=============================================================================
END MODULE Calibration
!=============================================================================

! $Log$
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
