! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Radiances ! Determine radiances for the GHz module
!=============================================================================

  USE MLSCommon, ONLY: r4, r8
  USE MLSL1Common, ONLY: GHzNum, FBchans, MBnum, MBchans, WFnum, WFchans, tau, &
       DACSnum, DACSchans, deflt_gain, deflt_zero, absZero_C, BandWidth, LO1
  USE MLSL1Utils, ONLY : GetIndexedAvg, Finite
  USE EngTbls, ONLY : Eng_MAF_T, CalTgtIndx, ReflecIndx, Reflec
  USE Calibration, ONLY : CalWin, limb_cnts, space_interp, target_interp, &
       space_err, target_err, Tsys, Cgain
  USE MLSL1Rad, ONLY : FBrad, MBrad, WFrad, DACSrad, RadPwr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CalcLimbRads

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CalcRadiance (limb_counts, space_counts, target_counts, &
       zero_counts, sum_w2, sum_wg2, space_P, target_P, baffle_P, bandwidth, &
       radNum, rad, rad_err, deflt_gain, use_deflt_gain, Tsys, gain, dacs)
!=============================================================================

    INTEGER :: radNum
    REAL(r8) :: limb_counts, space_counts, target_counts, zero_counts
    REAL(r8) :: sum_w2, sum_wg2
    REAL :: space_P, target_P, baffle_P, bandwidth, deflt_gain
    REAL :: rad, rad_err, gain, Tsys
    LOGICAL :: use_deflt_gain
    INTEGER, OPTIONAL :: dacs

    REAL :: baffle_L, baffle_S, baffle_T   ! Baffle radiances

    TYPE Eta_T
       CHARACTER(len=3) :: Name
       REAL(r4) :: Limb, Space, Target
    END TYPE Eta_T

    TYPE (Eta_T), PARAMETER :: Eta(0:4) = (/ &
         Eta_T ("R1A", 0.99598, 0.99587, 0.99575), &  ! Index '0'
         Eta_T ("R1B", 0.99344, 0.99317, 0.99274), &
         Eta_T ("R2 ", 0.99889, 0.99890, 0.99889), &
         Eta_T ("R3 ", 0.99929, 0.99928, 0.99929), &
         Eta_T ("R4 ", 0.99907, 0.99907, 0.99907)  /)

    rad = 0.0        ! nothing yet
    rad_err = -1.0   ! nothing yet

!! All baffle radiances same for now

    baffle_L = baffle_P
    baffle_S = baffle_P
    baffle_T = baffle_P

    IF (limb_counts == 0.0d0 .OR. space_counts == target_counts) RETURN

    IF (use_deflt_gain) THEN
       gain = deflt_gain
    ELSE
       gain =  (target_counts - space_counts) / &
            (Eta(radNum)%Target * target_P - Eta(radNum)%Space * space_P &
            - (1.0 - Eta(radNum)%Target) * baffle_T &
            + (1.0 - Eta(radNum)%Space) * baffle_S)
    ENDIF

    IF (gain == 0.0) RETURN  ! GAINS can be NEGATIVE!!!

! Calculate Tsys if it hasn't been done:

    IF (Tsys == 0.0) Tsys = (space_counts - zero_counts) / gain - space_P

    rad =  ((limb_counts - space_counts) / gain + Eta(radNum)%Space * space_P &
         - (1.0 - Eta(radNum)%Limb) * baffle_L + (1.0 - Eta(radNum)%Space) * &
         baffle_S) / Eta(radNum)%Limb

    IF (PRESENT (dacs)) THEN
       rad_err = (Tsys + rad)**2 + Tsys**2 * sum_w2 + rad**2 * &
            ((target_P + Tsys) / target_P)**2 + (1.0 + sum_w2) * sum_wg2
       rad_err = SQRT (rad_err / (bandwidth * tau))
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
       rad_err = SQRT (rad_err / (bandwidth * tau)) / gain
    ENDIF

  END SUBROUTINE CalcRadiance

!=============================================================================
  SUBROUTINE CalcNonLimbRad (Band, chan, RadNum, ReflecK, NonLimbRad)
!=============================================================================

    USE BandTbls, ONLY: SideBandFrac, SpilloverLoss, RadiometerLoss, BandFreq

    INTEGER, INTENT (IN) :: Band, chan, RadNum
    REAL, INTENT (IN) :: ReflecK
    REAL, INTENT (OUT) :: NonLimbRad

    REAL :: POA_L, POA_U, Spillover_L, Spillover_U
    REAL, PARAMETER :: GHzToHz = 1.0E09

    POA_L = RadPwr (BandFreq(Band)%Lower(chan)*GHzToHz, ReflecK)
    POA_U = RadPwr (BandFreq(Band)%Upper(chan)*GHzToHz, ReflecK)

    NonLimbRad = SideBandFrac(Band)%Lower(chan) * (1.0 - &
         RadiometerLoss(RadNum)%Ohmic) * POA_L
    NonLimbRad = NonLimbRad + SideBandFrac(Band)%Upper(chan) * (1.0 - &
         RadiometerLoss(RadNum)%Ohmic) * POA_U
    IF (Band < 32) THEN
       Spillover_L = SpilloverLoss(Band)%Lower(1,1)     ! Use only first value
       Spillover_U = SpilloverLoss(Band)%Upper(1,1)
    ELSE
       Spillover_L = SpilloverLoss(Band)%Lower(1,chan)  ! Use only first value
       Spillover_U = SpilloverLoss(Band)%Upper(1,chan)
    ENDIF
    Spillover_L = Spillover_L * RadiometerLoss(RadNum)%Spillover
    Spillover_U = Spillover_U * RadiometerLoss(RadNum)%Spillover

    NonLimbRad = NonLimbRad + SideBandFrac(Band)%Lower(chan) * &
         (1.0 - Spillover_L) * &
         RadiometerLoss(RadNum)%Ohmic * RadiometerLoss(RadNum)%Radiance
    NonLimbRad = NonLimbRad + SideBandFrac(Band)%Upper(chan) * &
         (1.0 - Spillover_U) * &
         RadiometerLoss(RadNum)%Ohmic * RadiometerLoss(RadNum)%Radiance

  END SUBROUTINE CalcNonLimbRad

!=============================================================================
  SUBROUTINE CalcLimbRads
!=============================================================================

    USE MLSL1Config, ONLY: MIFsGHz, L1Config
    USE MLSL1Common, ONLY: L1BFileInfo, Deflt_chi2
    USE MLSL1Utils, ONLY: Finite

    TYPE (Eng_MAF_T) :: EMAF

    CHARACTER(LEN=1) :: GHz_Cal_Type
    INTEGER :: bank, chan, MIF_index, MIF_index_MAX, radNum, bandNo, i
    INTEGER :: time_index, start_index, end_index, windex, CalWin_end
    REAL :: GHz_T1, GHz_T2, space_T, target_T, GHz_target_T, gain
    REAL :: ReflecAvg, NonLimbRad, Scf, Tcf
    REAL :: space_P, target_P, baffle_P ! Power per unit bandwidth
    REAL(r8) :: C_zero
    LOGICAL :: use_deflt_gains, do_chi2_err
    INTEGER, PARAMETER :: DACS_FB(4) = (/ 9, 2, 7, 1 /)

    use_deflt_gains = L1Config%Calib%UseDefaultGains

    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    CalWin_end = CalWin%MAFdata(CalWin%size)%end_index
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

    !    if (finite (GHz_T1)) print *, "GHzAmb avg: ", GHz_T1
    GHz_T2 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzCntl) - absZero_C
    !    if (finite (GHZ_T2)) print *, "GHzCntl avg: ", GHz_T2

    !  Temperature combinations for S and T from cf file:
!
    !    Note: no input T is same as "-"  and all final temps will be absolutes.
!
    !   Scf   Tcf   S   T
    !--------------------
    !    +     -    S   A
    !    +     0    S   C
    !    +     +    S   T
    !    -     -    A   S
    !    -     0    C   S
    !    0     0    A   C
    !    0     -    C   A

    Scf = L1Config%Calib%GHzSpaceTemp
    Tcf = L1Config%Calib%GHzTargetTemp

    IF (Scf > 0.0 .AND. Tcf < 0.0) THEN
       space_T = Scf
       GHz_target_T = GHz_T1
    ELSE IF (Scf > 0.0 .AND. NINT (Tcf) == 0) THEN
       space_T = Scf
       GHz_target_T = GHz_T2
    ELSE IF (Scf > 0.0 .AND. Tcf > 0.0) THEN
       space_T = Scf
       GHz_target_T = Tcf
    ELSE IF (Scf < 0.0 .AND. Tcf < 0.0) THEN
       space_T = GHz_T1
       GHz_target_T = ABS (Scf)
    ELSE IF (Scf < 0.0 .AND. NINT (Tcf) == 0) THEN
       space_T = GHz_T2
       GHz_target_T = ABS (Scf)
    ELSE IF (NINT (Scf) == 0 .AND. NINT (Tcf) == 0) THEN
       space_T = GHz_T1
       GHz_target_T = GHz_T2
    ELSE IF (NINT (Scf) == 0 .AND. Tcf < 0.0) THEN
       space_T = GHz_T2
       GHz_target_T = GHz_T1
    ENDIF

    ! Check if Temperatures are good

    IF (.NOT. Finite (space_T) .OR. .NOT. Finite (GHz_target_T) ) THEN
       DO bank = 1, GHzNum
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
    PRINT *, 'S/T temp: ', space_T, GHz_target_T
    Tsys%FB = 0.0
    Tsys%MB = 0.0
    Tsys%WF = 0.0
    Cgain%FB = 0.0
    Cgain%MB = 0.0
    Cgain%WF = 0.0
    MIF_index_MAX = MIFsGHz - 1

    ! Reflector temperatures:

    Reflec%Pri = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Pri) - absZero_C
    Reflec%Sec = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Sec) - absZero_C
    Reflec%Ter = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Ter) - absZero_C
    ReflecAvg = SUM((/ Reflec%Pri, Reflec%Sec, Reflec%Ter /)) / 3

    DO time_index = start_index, end_index   ! for every MIF in the MAF

       MIF_index = time_index - start_index  ! MIF # within the central MAF

       DO bank = 1, GHzNum

          ! Determine which Target temp to use!!

          target_T =  GHz_target_T
          BandNo = FBrad(bank)%bandno
          radNum = FBrad(bank)%signal%radiometerNumber
          IF (FBrad(bank)%signal%radiometerModifier == "A" .AND. &
               radNum == 1) radNum = 0   ! R1A
          space_P = radPwr (LO1(radNum), space_T)
          target_P = radPwr (LO1(radNum), target_T)
          baffle_P = radPwr (LO1(radNum), GHz_T1)
          do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

          IF (MIF_index <= MIF_index_MAX .AND. radNum < 5) THEN
             DO chan = 1, FBchans
                C_zero = deflt_zero%FB(chan,bank)
                CALL CalcRadiance (limb_cnts%FB(time_index,chan,bank), &
                     space_interp(MIF_index)%FB(chan,bank), &
                     target_interp(MIF_index)%FB(chan,bank), C_zero, &
                     space_err(MIF_index)%FB(chan,bank), &
                     target_err(MIF_index)%FB(chan,bank), &
                     space_P, target_P, baffle_P, BandWidth%FB(chan,bank), &
                     radNum, FBrad(bank)%value(chan,MIF_index+1), &
                     FBrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%FB(chan,bank), use_deflt_gains, &
                     Tsys%FB(chan,bank), gain)

                CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                     NonLimbRad)
                FBrad(bank)%value(chan,MIF_index+1) = &
                     FBrad(bank)%value(chan,MIF_index+1) - NonLimbRad

                IF (do_chi2_err) FBrad(bank)%precision(chan,MIF_index+1) = &
                     deflt_chi2%FB(chan,BandNo) * &
                     FBrad(bank)%precision(chan,MIF_index+1)

                IF (Cgain%FB(chan,bank) == 0.0 .AND. Tsys%FB(chan,bank) > 0.0) &
                     Cgain%FB(chan,bank) = gain

             ENDDO
          ENDIF

       ENDDO

       IF (MIF_index <= MIF_index_MAX) THEN
          DO bank = 1, MBnum

             bandNo = MBrad(bank)%bandno
             radNum = MBrad(bank)%signal%radiometerNumber
             space_P = radPwr (LO1(radNum), space_T)
             target_P = radPwr (LO1(radNum), target_T)
             baffle_P = radPwr (LO1(radNum), GHz_T1)
             do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

             DO chan = 1, MBchans
                C_zero = deflt_zero%MB(chan,bank)
                CALL CalcRadiance (limb_cnts%MB(time_index,chan,bank), &
                     space_interp(MIF_index)%MB(chan,bank), &
                     target_interp(MIF_index)%MB(chan,bank), C_zero, &
                     space_err(MIF_index)%MB(chan,bank), &
                     target_err(MIF_index)%MB(chan,bank), &
                     space_P, target_P, baffle_P, BandWidth%MB(chan,bank), &
                     radNum, MBrad(bank)%value(chan,MIF_index+1), &
                     MBrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%MB(chan,bank), use_deflt_gains, &
                     Tsys%MB(chan,bank), gain)

                CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                     NonLimbRad)
                MBrad(bank)%value(chan,MIF_index+1) = &
                     MBrad(bank)%value(chan,MIF_index+1) - NonLimbRad

                IF (do_chi2_err) MBrad(bank)%precision(chan,MIF_index+1) = &
                     deflt_chi2%MB(chan,bank) * &
                     MBrad(bank)%precision(chan,MIF_index+1)

                IF (Cgain%MB(chan,bank) == 0.0 .AND. Tsys%MB(chan,bank) > 0.0) &
                     Cgain%MB(chan,bank) = gain

             ENDDO
          ENDDO

          DO bank = 1, WFnum

             bandNo = WFrad(bank)%bandno
             radNum = WFrad(bank)%signal%radiometerNumber
             IF (WFrad(bank)%signal%radiometerModifier == "A" .AND. &
                  radNum == 1) radNum = 0   ! R1A
             space_P = radPwr (LO1(radNum), space_T)
             target_P = radPwr (LO1(radNum), target_T)
             baffle_P = radPwr (LO1(radNum), GHz_T1)
             do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

             DO chan = 1, WFchans
                C_zero = deflt_zero%WF(chan,bank)
                CALL CalcRadiance (limb_cnts%WF(time_index,chan,bank), &
                     space_interp(MIF_index)%WF(chan,bank), &
                     target_interp(MIF_index)%WF(chan,bank), C_zero, &
                     space_err(MIF_index)%WF(chan,bank), &
                     target_err(MIF_index)%WF(chan,bank), &
                     space_P, target_P, baffle_P, BandWidth%WF(chan,bank), &
                     radNum, WFrad(bank)%value(chan,MIF_index+1), &
                     WFrad(bank)%precision(chan,MIF_index+1), &
                     deflt_gain%WF(chan,bank), use_deflt_gains, &
                     Tsys%WF(chan,bank), gain)

                CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                     NonLimbRad)
                WFrad(bank)%value(chan,MIF_index+1) = &
                     WFrad(bank)%value(chan,MIF_index+1) - NonLimbRad

                IF (do_chi2_err) WFrad(bank)%precision(chan,MIF_index+1) = &
                     deflt_chi2%WF(chan,bank) * &
                     WFrad(bank)%precision(chan,MIF_index+1)

                IF (Cgain%WF(chan,bank) == 0.0 .AND. Tsys%WF(chan,bank) > 0.0) &
                     Cgain%WF(chan,bank) = gain

             ENDDO
          ENDDO

          IF (L1Config%Calib%CalibDACS) THEN

             DO bank = 1, DACSnum

                bandNo = DACSrad(bank)%bandno
                radNum = DACSrad(bank)%signal%radiometerNumber
                IF (DACSrad(bank)%signal%radiometerModifier == "A" .AND. &
                     radNum == 1) radNum = 0   ! R1A
                space_P = radPwr (LO1(radNum), space_T)
                target_P = radPwr (LO1(radNum), target_T)
                baffle_P = radPwr (LO1(radNum), GHz_T1)

                ! Tsys from the appropiate FB:

                IF (bank == 1 .AND. radNum == 1) THEN ! DACS 1 has been switched
                   Tsys%DACS(:,bank) = Tsys%FB(13,1)  ! Use FB 1's Tsys
                ELSE
                   Tsys%DACS(:,bank) = Tsys%FB(13,DACS_FB(bank))
                ENDIF

                DO chan = 1, DACSchans
                   C_zero = 0
                   CALL CalcRadiance (limb_cnts%DACS(time_index,chan,bank), &
                        space_interp(MIF_index)%DACS(chan,bank), &
                        target_interp(MIF_index)%DACS(chan,bank), C_zero, &
                        space_err(MIF_index)%DACS(chan,bank), &
                        target_err(MIF_index)%DACS(chan,bank), &
                        space_P, target_P, baffle_P, &
                        BandWidth%DACS(chan,bank), radNum, &
                        DACSrad(bank)%value(chan,MIF_index+1), &
                        DACSrad(bank)%precision(chan,MIF_index+1), &
                        deflt_gain%DACS(chan,bank), use_deflt_gains, &
                        Tsys%DACS(chan,bank), gain, dacs=1)

                   CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                        NonLimbRad)
                   DACSrad(bank)%value(chan,MIF_index+1) = &
                        DACSrad(bank)%value(chan,MIF_index+1) - NonLimbRad

                ENDDO
             ENDDO

          ENDIF

       ENDIF
    ENDDO

  END SUBROUTINE CalcLimbRads

!=============================================================================
END MODULE Radiances
!=============================================================================

! $Log$
! Revision 2.9  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.8  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.7  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.6  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.5  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.3  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.2  2001/09/10 16:17:56  perun
! Added CalMAFdata from Calibration module
!
! Revision 2.1  2001/02/23 20:55:04  perun
! Version 0.5 commit
!
