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
MODULE Radiances ! Determine radiances for the GHz module
!=============================================================================

  USE MLSCommon, ONLY: r4, r8
  USE MLSL1Common, ONLY: GHzNum, FBchans, MBnum, MBchans, WFnum, WFchans, tau, &
       DACSnum, DACSchans, deflt_gain, deflt_zero, absZero_C, BandWidth, LO1, &
       MaxMIFs, BandChanBad, limb_cnts, space_interp, target_interp, &
       target_err, slimb_interp, slimb_err, slimb_type, space_err
  USE MLSL1Utils, ONLY : GetIndexedAvg
  USE EngTbls, ONLY : Eng_MAF_T, CalTgtIndx, ReflecIndx, Reflec
  USE Calibration, ONLY : CalWin, Tsys, Cgain
  USE MLSL1Rad, ONLY : FBrad, MBrad, WFrad, DACSrad, RadPwr

  USE ALLOCATE_DEALLOCATE, ONLY : test_allocate, test_deallocate
  USE MLSSignalNomenclature, ONLY:GetFullMLSSignalName

  ! <whd> 
  !
  ! On the naming convention here. {FB,MB,WF,DACS}rad are pointers into an array
  ! named L1Brad of instances of the user-type MLSL1Rad::Radiance_T. You'll
  ! struggle in vain to find many references to this enclosing type. in
  ! particular, if you're looking for where the sub-field 'value' and
  ! 'precision' are set using the name 'L1Brad', you'll end up beating your head
  ! agains the wall. You won't find it. Instead, the individual portions are
  ! used in this routine, so you should look for FBrad(xx)%value, ..., etc.
  !
  ! Such practices makes maintenance harder, but there you have it. 
  ! Useful REs to find such things.
  !
  ! % grep -inH -P -e '\b(FB|MB|WF|DACS)Rad\(.*\)%value' *.f90
  !
  !</whd>

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CalcLimbRads

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CalcRadiance (limb_counts, space_counts, target_counts, &           ! in
       & zero_counts, sum_w2, sum_wg2, space_P, target_P, baffle_P, bandwidth, & ! in
       & radNum, bandNo, chanNo, &                                               ! in
       & rad, rad_err, &                                                         ! out
       & deflt_gain, use_deflt_gain, &                                           ! in
       & Tsys, gain, &                                                           ! out
       & slimb_type, slimb_counts, &                                             ! in
       & P_offset, &                                                             ! out
       & dacs)                                                                   ! in
!=============================================================================

    USE BandTbls, ONLY: GetEta_TSL

    INTEGER, INTENT (IN) :: bandNo, chanNo, radNum
    REAL(r8), INTENT(IN) :: limb_counts, space_counts, target_counts, zero_counts, &
         slimb_counts
    REAL(r8), INTENT(IN) :: sum_w2, sum_wg2
    REAL,INTENT(IN) :: space_P, target_P, baffle_P, bandwidth, deflt_gain
    REAL, INTENT(OUT) :: P_offset
    REAL,INTENT(OUT) :: rad, rad_err, gain, Tsys
    LOGICAL,INTENT(IN) :: use_deflt_gain, slimb_type
    INTEGER, INTENT(IN),OPTIONAL :: dacs

    ! Local variables

    REAL :: TsysD,baffle_L, baffle_S, baffle_T   ! Baffle radiances
    REAL :: eta_TSL(3)                     ! Target, Space and Limb eta values

    TYPE Rho_T
       CHARACTER(len=3) :: Name
       REAL(r4) :: Limb
    END TYPE Rho_T

    TYPE (Rho_T), PARAMETER :: Rho(0:4) = (/ &
         Rho_T ("R1A", 0.9938 * 0.9800), &  ! Index '0'
         Rho_T ("R1B", 0.9795 * 0.9753), &
         Rho_T ("R2 ", 0.9914 * 0.9945), &
         Rho_T ("R3 ", 0.9746 * 0.9953), &
         Rho_T ("R4 ", 0.9802 * 0.9888)  /)

    rad = 0.0        ! nothing yet
    rad_err = -1.0   ! nothing yet

!! All baffle radiances same for now

    baffle_L = baffle_P
    baffle_S = baffle_P
    baffle_T = baffle_P

    CALL GetEta_TSL (bandNo, chanNO, eta_TSL)

    IF (limb_counts == 0.0d0 .OR. space_counts == target_counts) RETURN

    IF (use_deflt_gain) THEN
       gain = deflt_gain
    ELSE
       ! Eq. 4.18 in ATB, p26
       gain =  (target_counts - space_counts) / &
            (eta_TSL(1) * target_P - eta_TSL(2) * space_P &
            - (1.0 - eta_TSL(1)) * baffle_T &
            + (1.0 - eta_TSL(2)) * baffle_S)
     ENDIF

    IF (gain == 0.0) RETURN  ! GAINS can be NEGATIVE!!!

! Calculate Tsys if it hasn't been done:

    ! See eq. 4.47 in the ATB, p 40.
    IF (Tsys == 0.0) Tsys = (space_counts - zero_counts) / gain - space_P

    IF (slimb_type) THEN
       rad = ((limb_counts - slimb_counts) / gain) / eta_TSL(3) + &
            Rho(radNum)%Limb * space_P
    ELSE
       !<whd> equation ATB 4.19, p 27</whd>
       rad = ((limb_counts - space_counts) / gain + eta_TSL(2)*space_P &
            - (1.0 - eta_TSL(3))*baffle_L + (1.0 - eta_TSL(2)) * &
            baffle_S) / eta_TSL(3)
    ENDIF

    !<whd> What is P_offset? Can't find any mention in the ATB </whd>
    P_offset = ((limb_counts - space_counts) / gain - space_P * &
         (eta_TSL(3) * Rho(radNum)%Limb - eta_TSL(2)) - &
         (1.0 - eta_TSL(3)) * baffle_L + (1.0 - eta_TSL(2)) * &
         baffle_S) / eta_TSL(2)

    IF (PRESENT (dacs)) THEN
       TsysD = Tsys / 0.87
       rad_err = (TsysD + rad)**2 + TsysD**2 * sum_w2 + rad**2 * &
            ((target_P + TsysD) / target_P)**2 + (1.0 + sum_w2) * sum_wg2
       rad_err = SQRT (rad_err / (bandwidth * tau))
    ELSE
       IF (slimb_type) THEN
          rad_err = &
               (limb_counts-zero_counts)**2 + (slimb_counts-zero_counts)**2 &
               * sum_w2 + (limb_counts-slimb_counts)**2 * &
               ((target_counts-zero_counts) / (target_counts-slimb_counts))**2 &
               * (1.0+sum_w2) * sum_wg2
       ELSE
          ! ATB, eq D.12. sum_w2,sum_wg2 is the sum-of-squares of the space and
          ! target interpolation weights, respectively
          rad_err = &
               (limb_counts-zero_counts)**2 + (space_counts-zero_counts)**2 &
               * sum_w2 + (limb_counts-space_counts)**2 * &
               ((target_counts-zero_counts) / (target_counts-space_counts))**2 &
               * (1.0+sum_w2) * sum_wg2
       ENDIF
       IF (rad_err < 0.0) THEN
          rad_err = -1.0
          RETURN
       ENDIF
       rad_err = SQRT (rad_err / (bandwidth * tau)) / gain
    ENDIF

  END SUBROUTINE CalcRadiance

!=============================================================================
  SUBROUTINE CalcNonLimbRad (Band, chan, RadNum, ReflecK, delrad, NonLimbRad)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE BandTbls, ONLY: SideBandFrac, SpilloverLoss, RadiometerLoss, BandFreq

    INTEGER, INTENT (IN) :: Band, chan, RadNum
    REAL, INTENT (IN) :: ReflecK, delrad
    REAL, INTENT (OUT) :: NonLimbRad


    ! Local Variables

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
         RadiometerLoss(RadNum)%Ohmic * (RadiometerLoss(RadNum)%Radiance + &
         delrad)
    NonLimbRad = NonLimbRad + SideBandFrac(Band)%Upper(chan) * &
         (1.0 - Spillover_U) * &
         RadiometerLoss(RadNum)%Ohmic * (RadiometerLoss(RadNum)%Radiance + &
         delrad)

    ! Scale based in user input in the cf file:

    NonLimbRad = NonLimbRad * L1Config%Calib%AntOffsetsScale

  END SUBROUTINE CalcNonLimbRad

!=============================================================================
  SUBROUTINE CalcLimbRads
!=============================================================================

    USE MLSL1Debug, ONLY: DebugControl, writeRadiancesInfo,SnoopRadiances
    USE SnoopMLSL1, ONLY: Snoop, L1SnoopOffering_T
    USE STRING_TABLE, ONLY: Create_String

    USE MLSL1Config, ONLY: MIFsGHz, L1Config
    USE MLSL1Common, ONLY: Deflt_chi2
    USE MLSL1Utils, ONLY: Finite
    USE SpectralBaseline, ONLY: CalcBaseline
    USE BandTbls, ONLY: GetDeltaRads, delrad_1_31, delrad_32_34

    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info

    ! MLSL1Debug and Snoop variables

    CHARACTER(len=256) SignalName
    CHARACTER(len=256) :: Location, Comment, msg
    REAL(r8),ALLOCATABLE :: TAI93(:) ! Time passed to snooper
    INTEGER status
    TYPE (L1SnoopOffering_T), Allocatable :: Offerings1(:),Offerings2(:)
    ! end MLSL1Debug


    TYPE (Eng_MAF_T) :: EMAF

    CHARACTER(LEN=1) :: GHz_Cal_Type

    INTEGER :: bank, chan, MIF_index, MIF_index_MAX, radNum, bandNo, Pnum
    INTEGER :: time_index, start_index, end_index, windex, CalWin_end
    REAL :: GHz_T1, GHz_T2, space_T, target_T, GHz_target_T, gain, rad_prec
    REAL :: ReflecAvg, NonLimbRad, Scf, Tcf, MIFprecSign(0:(MaxMIFs-1)), Pavg
    REAL :: space_P, target_P, baffle_P ! Power per unit bandwidth
    REAL(r8) :: C_zero
    LOGICAL :: use_deflt_gains, do_chi2_err, cntl_T
    INTEGER, PARAMETER :: DACS_FB(4) = (/ 9, 2, 7, 1 /)
    logical :: DeeBug



    print *,'In CalcLimbRads =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*'

    use_deflt_gains = L1Config%Calib%UseDefaultGains

    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    CalWin_end = CalWin%MAFdata(CalWin%size)%end_index
    EMAF = CalWin%MAFdata(windex)%EMAF
    MIFprecSign = CalWin%MAFdata(windex)%MIFprecSign


    WRITE ( msg,'("Central MAF no: ",i3)') EMAF%MAFno
    PRINT *,TRIM(msg)
    CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))


    CALL GetDeltaRads (CalWin%MAFdata(windex)%scGeodAngle)

    IF (ANY (INDEX (CalWin%MAFdata(windex)%scipkt(0:147)%GHz_sw_pos, "T") == &
         1)) THEN
      GHz_Cal_Type = "Primary"
    ELSE IF (ANY (INDEX (CalWin%MAFdata(windex)%scipkt(0:147)%GHz_sw_pos, "t") &
         == 1)) THEN
      GHz_Cal_Type = "Secondary"
    ELSE
      GHz_Cal_Type = "Unknown"
    ENDIF

    WRITE(msg,*) "GHz_Cal_Type = "//TRIM(GHz_Cal_Type)
    PRINT *,TRIM(msg)
    CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))

    GHz_T1 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzAmb) - absZero_C

    !    if (finite (GHz_T1)) print *, "GHzAmb avg: ", GHz_T1
    GHz_T2 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzCntl) - absZero_C &
         + 0.5   ! offset for R1(A/B), R2 and R3
    !    if (finite (GHZ_T2)) print *, "GHzCntl avg: ", GHz_T2
    !
    ! <whd> I have *no* idea what this table is about! </whd>
    !
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
    !
    !<whd> S=space, A=ambient?, C=control? And what does '0', '+' and '-' signify? </whd>

    Scf = L1Config%Calib%GHzSpaceTemp
    Tcf = L1Config%Calib%GHzTargetTemp
    cntl_T = .FALSE.

    ! <whd> Scf and Tcf are taken from the L1CF file and are used to select
    ! which temperatures to use in the following code. If Scf >0, use the value
    ! given by keyword `GHzSpaceTemp' from the CF file for the 'space'
    ! temperature, otherwise use the average of the GHz Ambient target that
    ! appears in the various field labeled 'GHzAmb...' in the engineering data
    ! stream. (that's what the calls to GetIndexedAverage above produce)
    ! </whd>

    IF (Scf > 0.0 .AND. Tcf < 0.0) THEN
      space_T = Scf
      GHz_target_T = GHz_T1
    ELSE IF (Scf > 0.0 .AND. NINT (Tcf) == 0) THEN
      space_T = Scf
      GHz_target_T = GHz_T2
      cntl_T = .TRUE.
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
      cntl_T = .TRUE.
    ELSE IF (NINT (Scf) == 0 .AND. Tcf < 0.0) THEN
      space_T = GHz_T2
      GHz_target_T = GHz_T1
    ENDIF

    ! Check if Temperatures are good. If either are non-finite, make them as bad
    ! (i.e. precision==-1)

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
      WRITE(msg,*) 'S/T temps are infinite! Nothing to do'
      PRINT *,TRIM(msg)
      Call MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
      RETURN
    ENDIF

    WRITE(msg,'("S/T temp: ", f5.2,"/",f7.3)') space_T, GHz_target_T
    PRINT *, trim(msg)
    CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))

    Tsys%FB = 0.0
    Tsys%MB = 0.0
    Tsys%WF = 0.0
    Cgain%FB = 0.0
    Cgain%MB = 0.0
    Cgain%WF = 0.0
    MIF_index_MAX = MIFsGHz - 1

    ! Reflector temperatures:

    Reflec%Pri = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Pri,DebugControl%Radiances) - absZero_C
    Reflec%Sec = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Sec,DebugControl%Radiances) - absZero_C
    Reflec%Ter = GetIndexedAvg (EMAF%eng%value, ReflecIndx%Ter,DebugControl%Radiances) - absZero_C

    ! Page 24 of ATB discusses taking this average
    ReflecAvg = SUM((/ Reflec%Pri, Reflec%Sec, Reflec%Ter /)) / 3

    ! Target temp with offset added

    target_T =  GHz_target_T

    CALL writeRadiancesInfo(GHzNum,      &
         &      FBchans,     &
         &      EMAF%MAFno,  &
         &      EMAF%TotalMAF, &
         &      start_index, &
         &      end_index,   &
         &      CalWin%central, &
         &      target_T,    &
         &      space_T,     &
         &      GHz_T1,      &
         &      GHz_T2)

    DO time_index = start_index, end_index   ! for every MIF in the MAF
      DeeBug = ( time_index == start_index )
      MIF_index = time_index - start_index  ! MIF # within the central MAF

      DO bank = 1, GHzNum
        ! CalWin%MAFdata(windex)%LimbAltFlag(MIF_index)%FB(1,bank)
        BandNo = FBrad(bank)%bandno
        radNum = FBrad(bank)%signal%radiometerNumber
        IF (FBrad(bank)%signal%radiometerModifier == "A" .AND. &
             radNum == 1) radNum = 0   ! R1A
        
        if ( DeeBug ) then
          CALL GetFullMLSSignalName(FBrad(bank)%signal, SignalName) ! Concatenate SD names
          print *, trim(SignalName)
          print *, 'Band ', BandNo
          print *, 'Radiometer ', RadNum
        endif
        ! IF (MIF_Index==0) THEN 
        !   CALL GetFullMLSSignalName(FBrad(bank)%signal, R2AName)
        !   PRINT *,"CalcLimbRads: FB: working on "//TRIM(R2AName)
        ! ENDIF
        space_P = radPwr (LO1(radNum), space_T)
        IF (cntl_T .AND. radNum == 4) THEN  ! Controlled Target
          target_P = radPwr (LO1(radNum), (target_T+0.1))  ! additional 0.1 K
        ELSE
          target_P = radPwr (LO1(radNum), target_T)
        ENDIF
        baffle_P = radPwr (LO1(radNum), GHz_T1)
        do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

        ! <whd> MIF_index_MAX = MIFsGHz-1 = 124 (data is 0 indexed, so there
        ! are 125 total). This is the number of MIFs in the output L1B
        ! data. As explained by Rick Cofield, this is the number of actual
        ! data (147-number of cal and otherwise unusable MIFs) There are 9 cal
        ! MIFs (4 space, 3 target) plus 1 MIF between and on either side =
        ! 10. Plus the .cf file specifies about 10 MIFs as 'discard'.
        !!
        ! Also, radNum==5 is the THz module and this code is for the GHz
        ! module radiances (I think) 
        !!
        !!
        ! </whd>


        ! <whd> I have to admit, though, this seems strange. There may be only
        ! 125 useful slots, but I don't see that they have to be the first
        ! 125! </whd>

        IF (MIF_index <= MIF_index_MAX .AND. radNum < 5) THEN
          if ( DEEBUG ) print*, 'Doing ', FBchans, ' FBchans'
          DO chan = 1, FBchans
            C_zero = deflt_zero%FB(chan,bank)
            ! <whd:question> 
            ! Why, in the following bit, is there code setting
            ! %value at MIF_index+1?
            ! </whd:question>
            CALL CalcRadiance (limb_cnts%FB(time_index,chan,bank), &   ! in
                 & space_interp(MIF_index)%FB(chan,bank), &            ! in
                 & target_interp(MIF_index)%FB(chan,bank), &           ! in
                 &  C_zero, &                                          ! in
                 & space_err(MIF_index)%FB(chan,bank), &               ! in
                 & target_err(MIF_index)%FB(chan,bank), &              ! in
                 & space_P, target_P, baffle_P, BandWidth%FB(chan,bank), & ! in
                 & radNum, BandNo, chan, &                                 ! in
                 & FBrad(bank)%value(chan,MIF_index+1), &                  ! out(why MIF_index+1?)
                 & rad_prec, &                                             ! out
                 & deflt_gain%FB(chan,bank), &                             ! in
                 & use_deflt_gains, &                                      ! in
                 & Tsys%FB(chan,bank), &                                   ! in/out (if == 0, 
                                        ! value is calculated)
                 & gain, &                                                 ! out (==deflt_gain 
                                        ! if use_deflt_gain)
                 & slimb_type%FB(chan,bank), &                             ! in
                 & slimb_interp(MIF_index)%FB(chan,bank), &                ! in
                 & FBrad(bank)%Poffset(chan,MIF_index+1))                  ! out


            ! <whd> If the band/channel is marked as `bad' in the
            ! L1CF,BandChanBad%Sign(Bandno, chan) will be < 0, so rad_prec
            ! will be < 0 and thus it will be marked as unusable.</whd>
            IF (rad_prec > 0.0) rad_prec = rad_prec * MIN ( &
                 MIFprecSign(MIF_index), BandChanBad%Sign(Bandno, chan))
            ! Set precisions < 0 when out-of-lock
            IF (rad_prec > 0.0 .and. RadNum == 2 .and. &
              & CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2) < 0. ) then
              rad_prec = rad_prec * CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2)
              if ( DEEBUG .and. chan == 1 ) print *, 'FBrad_prec ', rad_prec
            endif
            FBrad(bank)%precision(chan,MIF_index+1) = rad_prec

            CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                 delrad_1_31(BandNo), NonLimbRad)
            IF (MIF_index == 0) FBrad(bank)%ModelOffset(chan) = NonLimbRad

            IF (.NOT. slimb_type%FB(chan,bank)) THEN
              IF (Finite (NonLimbRad)) THEN
                FBrad(bank)%value(chan,MIF_index+1) = &
                     FBrad(bank)%value(chan,MIF_index+1) - NonLimbRad
              ELSE
                FBrad(bank)%value(chan,MIF_index+1) = 0.0
                FBrad(bank)%precision(chan,MIF_index+1) = -1.0
              ENDIF
            ENDIF

            IF (do_chi2_err) FBrad(bank)%precision(chan,MIF_index+1) = &
                 SQRT (deflt_chi2%FB(chan,BandNo)) * &
                 FBrad(bank)%precision(chan,MIF_index+1)

            IF (Cgain%FB(chan,bank) == 0.0 .AND. Tsys%FB(chan,bank) > 0.0) &
                 Cgain%FB(chan,bank) = gain

          ENDDO ! Loop over GHz banks
        ENDIF ! MIF is between 0 and MaxMIFs and not radiometer 5
      ENDDO ! loop over GHz channels

      IF (MIF_index <= MIF_index_MAX) THEN
          if ( DEEBUG ) print*, 'Doing ', MBNum, ' MBNum'
        DO bank = 1, MBnum
          ! IF (MIF_index == 0) THEN 
          !   CALL GetFullMLSSignalName(MBrad(bank)%signal, R2AName)
          !   PRINT *,"CalcLimbRads: MB: working on "//TRIM(R2AName)
          ! ENDIF
          bandNo = MBrad(bank)%bandno
          radNum = MBrad(bank)%signal%radiometerNumber
          space_P = radPwr (LO1(radNum), space_T)
          IF (cntl_T .AND. radNum == 4) THEN  ! Controlled Target
            target_P = radPwr (LO1(radNum), (target_T+0.1)) ! additional 0.1
          ELSE
            target_P = radPwr (LO1(radNum), target_T)
          ENDIF
          baffle_P = radPwr (LO1(radNum), GHz_T1)
          do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)

        if ( DeeBug ) then
          CALL GetFullMLSSignalName(MBrad(bank)%signal, SignalName) ! Concatenate SD names
          print *, trim(SignalName)
          print *, 'Band ', BandNo
          print *, 'Radiometer ', RadNum
        endif
          DO chan = 1, MBchans
            C_zero = deflt_zero%MB(chan,bank)
            CALL CalcRadiance (limb_cnts%MB(time_index,chan,bank), &
                 space_interp(MIF_index)%MB(chan,bank), &
                 target_interp(MIF_index)%MB(chan,bank), C_zero, &
                 space_err(MIF_index)%MB(chan,bank), &
                 target_err(MIF_index)%MB(chan,bank), &
                 space_P, target_P, baffle_P, BandWidth%MB(chan,bank), &
                 radNum, BandNo, chan,MBrad(bank)%value(chan,MIF_index+1), &
                 rad_prec, deflt_gain%MB(chan,bank), use_deflt_gains, &
                 Tsys%MB(chan,bank), gain, slimb_type%MB(chan,bank), &
                 slimb_interp(MIF_index)%MB(chan,bank), &
                 MBrad(bank)%Poffset(chan,MIF_index+1))

            IF (rad_prec > 0.0) rad_prec = rad_prec * MIN ( &
                 MIFprecSign(MIF_index),  BandChanBad%Sign(Bandno, chan))
            ! Set precisions < 0 when out-of-lock
            IF (rad_prec > 0.0 .and. RadNum == 2 .and. &
              & CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2) < 0. ) then
              rad_prec = rad_prec * CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2)
              if ( DEEBUG .and. chan == 1 ) print *, 'MBrad_prec ', rad_prec
            endif
            MBrad(bank)%precision(chan,MIF_index+1) = rad_prec

            CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                 delrad_1_31(BandNo), NonLimbRad)
            IF (MIF_index == 0) MBrad(bank)%ModelOffset(chan) = NonLimbRad

            IF (.NOT. slimb_type%MB(chan,bank)) THEN
              IF (Finite (NonLimbRad)) THEN
                MBrad(bank)%value(chan,MIF_index+1) = &
                     MBrad(bank)%value(chan,MIF_index+1) - NonLimbRad
              ELSE
                MBrad(bank)%value(chan,MIF_index+1) = 0.0
                MBrad(bank)%precision(chan,MIF_index+1) = -1.0
              ENDIF
            ENDIF

            IF (do_chi2_err) MBrad(bank)%precision(chan,MIF_index+1) = &
                 SQRT (deflt_chi2%MB(chan,bank)) * &
                 MBrad(bank)%precision(chan,MIF_index+1)

            IF (Cgain%MB(chan,bank) == 0.0 .AND. Tsys%MB(chan,bank) > 0.0) &
                 Cgain%MB(chan,bank) = gain

          ENDDO
        ENDDO

        if ( DEEBUG ) print*, 'Doing ', WFNum, ' WFNum'
        DO bank = 1, WFnum
          bandNo = WFrad(bank)%bandno
          ! IF (MIF_index == 0 ) THEN 
          !   CALL GetFullMLSSignalName(WFrad(bank)%signal, R2AName)
          !   PRINT *,"CalcLimbRads: WF: working on "//TRIM(R2AName)
          ! ENDIF
          radNum = WFrad(bank)%signal%radiometerNumber
          IF (WFrad(bank)%signal%radiometerModifier == "A" .AND. &
               radNum == 1) radNum = 0   ! R1A
          space_P = radPwr (LO1(radNum), space_T)
          target_P = radPwr (LO1(radNum), target_T)
          baffle_P = radPwr (LO1(radNum), GHz_T1)
          do_chi2_err = L1Config%Output%EnableChi2Err(BandNo)
        
        if ( DeeBug ) then
          CALL GetFullMLSSignalName(WFrad(bank)%signal, SignalName) ! Concatenate SD names
          print *, trim(SignalName)
          print *, 'Band ', BandNo
          print *, 'Radiometer ', RadNum
        endif

          DO chan = 1, WFchans
            C_zero = deflt_zero%WF(chan,bank)
            CALL CalcRadiance (limb_cnts%WF(time_index,chan,bank), &
                 space_interp(MIF_index)%WF(chan,bank), &
                 target_interp(MIF_index)%WF(chan,bank), C_zero, &
                 space_err(MIF_index)%WF(chan,bank), &
                 target_err(MIF_index)%WF(chan,bank), &
                 space_P, target_P, baffle_P, BandWidth%WF(chan,bank), &
                 radNum, BandNo, chan,WFrad(bank)%value(chan,MIF_index+1), &
                 rad_prec, deflt_gain%WF(chan,bank), use_deflt_gains, &
                 Tsys%WF(chan,bank), gain, slimb_type%WF(chan,bank), &
                 slimb_interp(MIF_index)%WF(chan,bank), &
                 WFrad(bank)%Poffset(chan,MIF_index+1))

            IF (rad_prec > 0.0) rad_prec = rad_prec * MIN ( &
                 MIFprecSign(MIF_index), BandChanBad%Sign(Bandno, chan))
            ! Set precisions < 0 when out-of-lock
            IF (rad_prec > 0.0 .and. RadNum == 2 .and. &
              & CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2) < 0. ) then
              rad_prec = rad_prec * CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2)
              if ( DEEBUG .and. chan == 1 ) print *, 'WFrad_prec ', rad_prec
            endif
            WFrad(bank)%precision(chan,MIF_index+1) = rad_prec

            CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                 delrad_32_34(BandNo,chan), NonLimbRad)
            IF (MIF_index == 0) WFrad(bank)%ModelOffset(chan) = NonLimbRad

            IF (.NOT. slimb_type%WF(chan,bank)) THEN
              IF (Finite (NonLimbRad)) THEN
                WFrad(bank)%value(chan,MIF_index+1) = &
                     WFrad(bank)%value(chan,MIF_index+1) - NonLimbRad
              ELSE
                WFrad(bank)%value(chan,MIF_index+1) = 0.0
                WFrad(bank)%precision(chan,MIF_index+1) = -1.0
              ENDIF
            ENDIF

            IF (do_chi2_err) WFrad(bank)%precision(chan,MIF_index+1) = &
                 SQRT (deflt_chi2%WF(chan,bank)) * &
                 WFrad(bank)%precision(chan,MIF_index+1)

            IF (Cgain%WF(chan,bank) == 0.0 .AND. Tsys%WF(chan,bank) > 0.0) &
                 Cgain%WF(chan,bank) = gain

          ENDDO
        ENDDO

        IF (L1Config%Calib%CalibDACS) THEN

          DO bank = 1, DACSnum
            ! IF (MIF_index == 0) THEN 
            !   CALL GetFullMLSSignalName(DACSrad(bank)%signal, R2AName)
            !   PRINT *,"CalcLimbRads: DACS: working on "//TRIM(R2AName)
            ! ENDIF
            bandNo = DACSrad(bank)%bandno
            radNum = DACSrad(bank)%signal%radiometerNumber
            IF (DACSrad(bank)%signal%radiometerModifier == "A" .AND. &
                 radNum == 1) radNum = 0   ! R1A
            space_P = radPwr (LO1(radNum), space_T)
            target_P = radPwr (LO1(radNum), target_T)
            baffle_P = radPwr (LO1(radNum), GHz_T1)

        if ( DeeBug ) then
          CALL GetFullMLSSignalName(DACSrad(bank)%signal, SignalName) ! Concatenate SD names
          print *, trim(SignalName)
          print *, 'Band ', BandNo
          print *, 'Radiometer ', RadNum
        endif
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
                   BandWidth%DACS(chan,bank), radNum, BandNo, chan, &
                   DACSrad(bank)%value(chan,MIF_index+1), &
                   rad_prec, deflt_gain%DACS(chan,bank), use_deflt_gains, &
                   Tsys%DACS(chan,bank), gain, slimb_type%DACS(chan,bank),&
                   slimb_interp(MIF_index)%DACS(chan,bank), &
                   DACSrad(bank)%Poffset(chan,MIF_index+1), dacs=1)

              IF (rad_prec > 0.0) rad_prec = rad_prec * MIN ( &
                   MIFprecSign(MIF_index), BandChanBad%Sign(Bandno, chan))
              ! Set precisions < 0 when out-of-lock
              IF (rad_prec > 0.0 .and. RadNum == 2 .and. &
                & CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2) < 0. ) then
                rad_prec = rad_prec * CalWin%MAFdata(windex)%RnPrecSign(MIF_index, 2)
                if ( DEEBUG .and. chan == 1 ) print *, 'DACSrad_prec ', rad_prec
              endif
                DACSrad(bank)%precision(chan,MIF_index+1) = rad_prec

              CALL CalcNonLimbRad (BandNo, chan, radNum, ReflecAvg, &
                   delrad_1_31(BandNo), NonLimbRad)
              IF (MIF_index == 0) DACSrad(bank)%ModelOffset(chan) = &
                   NonLimbRad

              IF (.NOT. slimb_type%DACS(chan,bank)) THEN
                IF (Finite (NonLimbRad)) THEN
                  DACSrad(bank)%value(chan,MIF_index+1) = &
                       DACSrad(bank)%value(chan,MIF_index+1) - NonLimbRad
                ELSE
                  DACSrad(bank)%value(chan,MIF_index+1) = 0.0
                  DACSrad(bank)%precision(chan,MIF_index+1) = -1.0
                ENDIF
              ENDIF

            ENDDO
          ENDDO
        ENDIF ! if doing DACS
      ENDIF ! is MIF < MIF_Max
    ENDDO ! loop over MIFs

    IF (SnoopRadiances) THEN
      ! <<<start here >>>
      Location="Radiances::CalcLimRads"
      Comment="Just after FB radiances calculation"
      PRINT *,TRIM(Location)//' : '//TRIM(Comment)

      ALLOCATE(Offerings1(GHzNum), stat=status)

      DO bank=1,GHzNum
        CALL GetFullMLSSignalName(FBrad(bank)%signal, SignalName)
        Offerings1(bank)%name=SignalName
        Offerings1(bank)%rank=2
        Offerings1(bank)%dimensions=(/FBChans,MIFsGHz,0/)
        Offerings1(bank)%R2Value => FBrad(bank)%value
      END DO

      CALL Snoop( TRIM(location), Offerings1, Offerings2, TRIM(comment))
    ENDIF
    ! Avg Poffsets:

    DO bank = 1, GHzNum
      DO chan = 1, FBchans
        Pavg = 0.0
        Pnum = 0
        DO MIF_index = 0, min( (MaxMIFs-1), MIFsGHz-1 )
          IF (CalWin%MAFdata(windex)%LimbAltFlag(MIF_index)%FB(chan,bank) &
               .AND. FBrad(bank)%precision(chan,MIF_index+1) > 0.0) THEN
            Pnum = Pnum + 1
            Pavg = Pavg + FBrad(bank)%Poffset(chan,MIF_index+1)
          ENDIF
        ENDDO
        IF (Pnum > 0) Pavg = Pavg / Pnum

        ! Save avg in first MIF to output in DIAG file:

        FBrad(bank)%Poffset(chan,1) = Pavg

      ENDDO
    ENDDO

    DO bank = 1, MBnum
      DO chan = 1, MBchans
        Pavg = 0.0
        Pnum = 0
        DO MIF_index = 0, min( (MaxMIFs-1), MIFsGHz-1 )
          IF (CalWin%MAFdata(windex)%LimbAltFlag(MIF_index)%MB(chan,bank) &
               .AND. MBrad(bank)%precision(chan,MIF_index+1) > 0.0) THEN
            Pnum = Pnum + 1
            Pavg = Pavg + MBrad(bank)%Poffset(chan,MIF_index+1)
          ENDIF
        ENDDO
        IF (Pnum > 0) Pavg = Pavg / Pnum

        ! Save avg in first MIF to output in DIAG file:

        MBrad(bank)%Poffset(chan,1) = Pavg

      ENDDO
    ENDDO

    DO bank = 1, WFnum
      DO chan = 1, WFchans
        Pavg = 0.0
        Pnum = 0
        DO MIF_index = 0, min( (MaxMIFs-1), MIFsGHz-1 )
          IF (CalWin%MAFdata(windex)%LimbAltFlag(MIF_index)%WF(chan,bank) &
               .AND. WFrad(bank)%precision(chan,MIF_index+1) > 0.0) THEN
            Pnum = Pnum + 1
            Pavg = Pavg + WFrad(bank)%Poffset(chan,MIF_index+1)
          ENDIF
        ENDDO
        IF (Pnum > 0) Pavg = Pavg / Pnum

        ! Save avg in first MIF to output in DIAG file:

        WFrad(bank)%Poffset(chan,1) = Pavg

      ENDDO
    ENDDO

    DO bank = 1, DACSnum
      DO chan = 1, DACSchans
        Pavg = 0.0
        Pnum = 0
        DO MIF_index = 0, min( (MaxMIFs-1), MIFsGHz-1 )
          IF (CalWin%MAFdata(windex)%LimbAltFlag(MIF_index)%DACS(chan,bank) &
               .AND. DACSrad(bank)%precision(chan,MIF_index+1) > 0.0) THEN
            Pnum = Pnum + 1
            Pavg = Pavg + DACSrad(bank)%Poffset(chan,MIF_index+1)
          ENDIF
        ENDDO
        IF (Pnum > 0) Pavg = Pavg / Pnum

        ! Save avg in first MIF to output in DIAG file:

        DACSrad(bank)%Poffset(chan,1) = Pavg

      ENDDO
    ENDDO

    ! Spectral Baseline:

    CALL CalcBaseline

  END SUBROUTINE CalcLimbRads


!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE Radiances
!=============================================================================

! $Log$
! Revision 2.24  2024/10/10 20:19:35  pwagner
! Set R2 precisions negative when out-of-lock
!
! Revision 2.23  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.22.4.2  2016/03/14 19:51:24  whdaffer
! Most of the work is to eliminate cicular references between MLSL1Debug
! and Calibration. To resolve this, I've moved the inclusion of
! Calibration.f9h and the definition of some 10 variables from
! Calibration.f90 to MLSL1Common.f90. Radiances, MLSL1Debug and
! Calibration will get those types, variable from MLSL1Common. Also, use
! machines.f90 to get the definition of usleep used in SnoopMLSL1
!
! Revision 2.22.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.22  2015/01/13 18:42:17  pwagner
! Avoid blowing past upper bounds on precision arrays
!
! Revision 2.21  2007/02/09 15:05:59  perun
! Always calculate P_offset
!
! Revision 2.20  2006/09/28 16:17:06  perun
! Save only one ModelOffset per MAF and calculate average slimb view Poffsets
!
! Revision 2.19  2006/08/18 15:53:19  perun
! Replace slimb_err with space_err to correct rad_err calculation
!
! Revision 2.18  2006/06/14 13:48:44  perun
!  Adjust radiances using the delrad values from the BandTbls
!
! Revision 2.17  2006/04/05 18:09:32  perun
! Remove unused variables
!
! Revision 2.16  2006/03/24 15:17:06  perun
! Add calculation of radiances using the Space/Limb view interpolations
!
! Revision 2.15  2005/12/09 16:39:58  perun
! Update port baffle transmission values
!
! Revision 2.14  2005/08/24 15:52:54  perun
! Set rads and precs to fill values when nonlimb radiances are unavailable
!
! Revision 2.13  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.12  2005/05/02 16:06:38  perun
! Use AntOffsetsScale for scaling NonLimbRad
!
! Revision 2.11  2004/11/15 16:50:06  perun
! Adjust controlled target temperature per RFJ
!
! Revision 2.10  2004/11/10 15:40:40  perun
! Adjust precision based on flag; change DACS precision method; call baseline
! calculation
!
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
