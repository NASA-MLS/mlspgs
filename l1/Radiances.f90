! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Radiances ! Determine radiances
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
       DACSnum, DACSchans, Chan_R_T, deflt_gain, deflt_zero, MaxMIFs, &
       absZero_C, BandWidth, tau, LO1
  USE MLSL1Utils, ONLY : GetIndexedAvg, Finite
  USE EngTbls, ONLY : Eng_MAF_T, CalTgtIndx
  USE Calibration, ONLY : CalWin, limb_time, limb_counts, space_interp, &
       target_interp, space_err, target_err, space_weight, target_weight
  USE MLSL1Rad, ONLY : FBrad, MBrad, WFrad, DACSrad

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (Chan_R_T) :: Tsys, Cgain

CONTAINS

  FUNCTION RadPwr (Hz, T) RESULT (P)

    USE MLSL1Common, ONLY: boltz, planck

!! Calculate radiant power per unit bandwidth

    REAL :: Hz    !! frequency in Hz
    REAL :: T     !! temperature in Kelvin
    REAL :: P     !! radiant power

    P = (planck * Hz) / (boltz * (exp ((planck * Hz) / (boltz * T)) - 1.0))

  END FUNCTION RadPwr

  SUBROUTINE CalcRadiance (limb_counts, space_counts, target_counts, &
       zero_counts, sum_w2, sum_wg2, space_P, target_P, bandwidth, &
       rad, rad_err, deflt_gain, use_deflt_gain, Tsys, gain, dacs)

    REAL(r8) :: limb_counts, space_counts, target_counts, zero_counts
    REAL(r8) :: sum_w2, sum_wg2
    REAL :: space_P, target_P, bandwidth, rad, rad_err, deflt_gain, Tsys
    LOGICAL :: use_deflt_gain
    INTEGER, OPTIONAL :: dacs

    REAL :: gain

    rad = 0.0        ! nothing yet
    rad_err = -1.0   ! nothing yet

    IF (limb_counts == 0.0d0 .OR. space_counts == target_counts) RETURN

    IF (use_deflt_gain) THEN
       gain = deflt_gain
    ELSE
       gain =  (target_counts - space_counts) / (target_P - space_P)
    ENDIF

    IF (gain == 0.0) RETURN  ! GAINS can be NEGATIVE!!!

! Calculate Tsys if it hasn't been done:

    IF (Tsys == 0.0) Tsys = (space_counts - zero_counts) / gain - space_P

    rad =  (limb_counts - space_counts) / gain + space_P

    IF (PRESENT (dacs)) THEN
       rad_err = (Tsys + rad)**2 + Tsys**2 * sum_w2 + rad**2 * &
            ((target_P + Tsys) / target_P)**2 + (1.0 + sum_w2) * sum_wg2
       rad_err = sqrt (rad_err / (bandwidth * tau))

    ELSE
       rad_err = &
            (limb_counts - zero_counts)**2 + (space_counts - zero_counts)**2 &
            * sum_w2 + (limb_counts - space_counts)**2 * &
            ((target_counts - zero_counts) / (target_counts-space_counts))**2 &
            * (1.0 + sum_w2) * sum_wg2
       IF (rad_err < 0.0) THEN
          rad_err = -1.0
          RETURN
       ENDIF
       rad_err = sqrt (rad_err / (bandwidth * tau)) / gain
    ENDIF

  END SUBROUTINE CalcRadiance

  SUBROUTINE CalcLimbRads

    USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz, L1Config
    USE MLSL1Common, ONLY: L1BFileInfo

    TYPE (Eng_MAF_T) :: EMAF
    TYPE ValPrec_T
       REAL :: value
       REAL :: precision
    END TYPE ValPrec_T
    TYPE Rad_T
       TYPE (ValPrec_T) :: FB(FBchans,FBnum)          ! standard filter banks
       TYPE (ValPrec_T) :: MB(MBchans,MBnum)          ! mid-band filter banks
       TYPE (ValPrec_T) :: WF(WFchans,WFnum)          ! wide filters
       TYPE (ValPrec_T) :: DACS(DACSchans,DACSnum)    ! DACS filters
    END TYPE Rad_T
    TYPE (Rad_T) :: Rad(0:MaxMIFs-1)

    CHARACTER(LEN=1) :: GHz_Cal_Type
    INTEGER :: i, j, bank, chan, MIF_index, MIF_index_MAX, radNum
    INTEGER :: time_index, start_index, end_index, windex, CalWin_end
    REAL :: GHz_T1, GHz_T2, THz_T1, space_T, target_T, GHz_target_T, gain, temp
    REAL :: space_P, target_P ! Power per unit bandwidth
    REAL(r8) :: C_zero
    LOGICAL :: use_deflt_gains
    INTEGER, PARAMETER :: DACS_FB(4) = (/ 9, 2, 7, 1 /) ! TEMP!!!

    use_deflt_gains = L1Config%Calib%UseDefaultGains

    windex = CalWin%central
    EMAF = CalWin%MAFdata(windex)%EMAF

    IF (ANY (INDEX (CalWin%MAFdata(windex)%scipkt(0:147)%GHz_sw_pos, "T") == &
         1)) THEN
       GHz_Cal_Type = "Primary"
    ELSE IF (ANY (INDEX (CalWin%MAFdata(windex)%scipkt(0:147)%GHz_sw_pos, "t") &
         == 1)) THEN
       GHz_Cal_Type = "Secondary"
    ELSE
       GHz_Cal_Type = "Unknown"
    ENDIF

    GHz_T1 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzAmb) - absZero_C
!    print *, "GHzAmb avg: ", GHz_T1
    GHz_T2 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzCntl) - absZero_C
!    print *, "GHzCntl avg: ", GHz_T2
    THz_T1 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%THzAmb) - absZero_C
    if (finite(thz_t1)) print *, "THzAmb avg: ", THz_T1

    IF (L1Config%Calib%TargetTemp > 0.0) THEN  ! Use the config input
       GHz_target_T = L1Config%Calib%TargetTemp
    ELSE
       IF (GHz_Cal_Type == "P") THEN
          GHz_target_T = GHz_T1
       ELSE IF (GHz_Cal_Type == "S") THEN
          GHz_target_T = GHz_T2
       ELSE
          print *, "NO CAL TARGET TEMP!"
          stop   !what else can I do?
       ENDIF
    ENDIF

    space_T = L1Config%Calib%SpaceTemp

!! Special test config

    IF (space_T < 0.0 .AND. L1Config%Calib%TargetTemp <= 0.0) THEN
       GHz_target_T = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzCntl) - &
            absZero_C
       space_T = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzAmb) - absZero_C
    ENDIF

!! Case for linearity (switch T and S):

    IF (space_T < 0.0) THEN
       temp = GHz_target_T
       GHz_target_T = -space_T
       space_T = temp
    ENDIF

! Check if Temperatures are good

    IF (.NOT. Finite (space_T) .OR. .NOT. Finite (GHz_target_T) ) THEN
       DO bank = 1, FBnum
          FBrad(bank)%value = 0.0
          FBrad(bank)%precision = -1.0
       ENDDO
       DO bank = 1, MBnum
          MBrad(bank)%value = 0.0
          MBrad(bank)%precision = -1.0
       ENDDO
       DO bank = 1, WFnum
          WFrad(bank)%value = 0.0
          WFrad(bank)%precision = -1.0
       ENDDO
       DO bank = 1, DACSnum
          DACSrad(bank)%value = 0.0
          DACSrad(bank)%precision = -1.0
       ENDDO
       RETURN
    ENDIF
!print *, 'S/T temp: ', space_T, GHz_target_T
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    CalWin_end = CalWin%MAFdata(CalWin%size)%end_index
    Tsys%FB = 0.0
    Tsys%MB = 0.0
    Tsys%WF = 0.0
    Cgain%FB = 0.0
    Cgain%MB = 0.0
    Cgain%WF = 0.0

    DO time_index = start_index, end_index   ! for every MIF in the MAF

       MIF_index = time_index - start_index  ! MIF # within the central MAF

       DO bank = 1, FBnum

! Determine which Target temp to use!!

          target_T =  GHz_target_T

          radNum = FBrad(bank)%signal%radiometerNumber

          space_P = radPwr (LO1(radNum), space_T)
          target_P = radPwr (LO1(radNum), target_T)

! Determine MAX allowable MIF index

          IF (bank < 15) THEN
             MIF_index_MAX = MIFsGHz - 1
          ELSE
             MIF_index_MAX = MIFsTHz - 1
          ENDIF

          IF (MIF_index <= MIF_index_MAX) THEN
             DO chan = 1, FBchans
                C_zero = deflt_zero%FB(chan,bank)
                CALL CalcRadiance (limb_counts(time_index)%FB(chan,bank), &
                     space_interp(MIF_index)%FB(chan,bank), &
                     target_interp(MIF_index)%FB(chan,bank), C_zero, &
                     space_err(MIF_index)%FB(chan,bank), &
                     target_err(MIF_index)%FB(chan,bank), &
                     space_P, target_P, BandWidth%FB(chan,bank), &
                     FBrad(bank)%value(chan,MIF_index+1), &
                     FBrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%FB(chan,bank), use_deflt_gains, &
                     Tsys%FB(chan,bank), gain)
                IF (Cgain%FB(chan,bank) == 0.0 .AND. Tsys%FB(chan,bank) > 0.0) &
                     Cgain%FB(chan,bank) = gain

             ENDDO
          ENDIF

       ENDDO

       MIF_index_MAX = MIFsGHz - 1

       IF (MIF_index <= MIF_index_MAX) THEN
          DO bank = 1, MBnum

             radNum = MBrad(bank)%signal%radiometerNumber
             space_P = radPwr (LO1(radNum), space_T)
             target_P = radPwr (LO1(radNum), target_T)

             DO chan = 1, MBchans
                C_zero = deflt_zero%MB(chan,bank)
                CALL CalcRadiance (limb_counts(time_index)%MB(chan,bank), &
                     space_interp(MIF_index)%MB(chan,bank), &
                     target_interp(MIF_index)%MB(chan,bank), C_zero, &
                     space_err(MIF_index)%MB(chan,bank), &
                     target_err(MIF_index)%MB(chan,bank), &
                     space_P, target_P, BandWidth%MB(chan,bank), &
                     MBrad(bank)%value(chan,MIF_index+1), &
                     MBrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%MB(chan,bank), use_deflt_gains, &
                     Tsys%MB(chan,bank), gain)
                IF (Cgain%MB(chan,bank) == 0.0 .AND. Tsys%MB(chan,bank) > 0.0) &
                     Cgain%MB(chan,bank) = gain
                
             ENDDO
          ENDDO

          DO bank = 1, WFnum

             radNum = WFrad(bank)%signal%radiometerNumber
             space_P = radPwr (LO1(radNum), space_T)
             target_P = radPwr (LO1(radNum), target_T)

             DO chan = 1, WFchans
                C_zero = deflt_zero%WF(chan,bank)
                CALL CalcRadiance (limb_counts(time_index)%WF(chan,bank), &
                     space_interp(MIF_index)%WF(chan,bank), &
                     target_interp(MIF_index)%WF(chan,bank), C_zero, &
                     space_err(MIF_index)%WF(chan,bank), &
                     target_err(MIF_index)%WF(chan,bank), &
                     space_P, target_P, BandWidth%WF(chan,bank), &
                     WFrad(bank)%value(chan,MIF_index+1), &
                     WFrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%WF(chan,bank), use_deflt_gains, &
                     Tsys%WF(chan,bank), gain)
                IF (Cgain%WF(chan,bank) == 0.0 .AND. Tsys%WF(chan,bank) > 0.0) &
                     Cgain%WF(chan,bank) = gain
                
             ENDDO
          ENDDO

          IF (L1Config%Calib%CalibDACS) THEN

             BandWidth%DACS = 0.15 * 1.0e06    ! Bandwidths for all DACS
             DO bank = 1, DACSnum

                radNum = DACSrad(bank)%signal%radiometerNumber
                space_P = radPwr (LO1(radNum), space_T)
                target_P = radPwr (LO1(radNum), target_T)
                Tsys%DACS(:,bank) = Tsys%FB(13,DACS_FB(bank))   !TEMP !!!

                DO chan = 1, DACSchans
                   C_zero = 0
                   CALL CalcRadiance (limb_counts(time_index)%DACS(chan,bank), &
                        space_interp(MIF_index)%DACS(chan,bank), &
                        target_interp(MIF_index)%DACS(chan,bank), C_zero, &
                        space_err(MIF_index)%DACS(chan,bank), &
                        target_err(MIF_index)%DACS(chan,bank), &
                        space_P, target_P, BandWidth%DACS(chan,bank), &
                        DACSrad(bank)%value(chan,MIF_index+1), &
                        DACSrad(bank)%precision(chan,MIF_index+1), &
                        deflt_gain%DACS(chan,bank), use_deflt_gains, &
                        Tsys%DACS(chan,bank), gain, dacs=1)

                ENDDO
             ENDDO

          ENDIF

       ENDIF
    ENDDO

! Write Diags

    WRITE (L1BFileInfo%DiagId) CalWin%MAFdata(windex)%SciPkt(0)%MAFno
    WRITE (L1BFileInfo%DiagId) Tsys%FB, Tsys%MB, Tsys%WF
    WRITE (L1BFileInfo%DiagId) Cgain%FB, Cgain%MB, Cgain%WF

  END SUBROUTINE CalcLimbRads

!=============================================================================
END MODULE Radiances
!=============================================================================

! $Log$
! Revision 2.3  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.2  2001/09/10 16:17:56  perun
! Added CalMAFdata from Calibration module
!
! Revision 2.1  2001/02/23 20:55:04  perun
! Version 0.5 commit
!
