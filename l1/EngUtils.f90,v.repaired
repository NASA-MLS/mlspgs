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
MODULE EngUtils   ! Engineering utilities
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Utils, ONLY: QNan, BigEndianStr, Finite
  use MLSFillValues, only: isNaN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: NextEngMAF

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

CONTAINS

  FUNCTION Cal_freq (meas_freq, hi_freq, lo_freq, hi_cal, lo_cal) RESULT (freq)

    !! Calibrate the measured frequency against the calibration frequencies and
    !! their corresponding parameter values.  

    !! See Eq. 5.3, sec 5 in the _MLS Level Algorithmic Theoretical Basis_
    !! (ATB) version 2.0)

    !--------Arguments--------!
    INTEGER, INTENT(IN) :: meas_freq  ! measured frequency
    INTEGER, INTENT(IN) :: hi_freq    ! "high" calibration frequency
    INTEGER, INTENT(IN) :: lo_freq    ! "low" calibration frequency
    REAL, INTENT(IN) :: hi_cal        ! "high" calibration parameter
    REAL, INTENT(IN) :: lo_cal        ! "low" calibration parameter

    REAL :: freq                      ! calculated frequency

    IF ((lo_freq == 0) .OR. (hi_freq == 0) .OR. (meas_freq == 0.0) .OR. &
         (hi_freq == lo_freq)) THEN
       freq = QNan()
    ELSE
       freq = (REAL(meas_freq - lo_freq) / (hi_freq - lo_freq)) * &
            & (hi_cal - lo_cal) + lo_cal
    ENDIF

  END FUNCTION Cal_freq

  FUNCTION Therm_temp (rin) RESULT (T)

    !! See Eq.s 5.5, 5.6, Sec 5.1.3, page 49 of the MLS Level 1 ATB document
    !! (rev 2.0)
    !!

    !! Convert thermistor (ohms) to temperature (C) */
    !!
    !! The thermistor is a YSI 44906 in parallel with a 4.99 K ohm resistor.
    !!
    !! input:
    !!         rin: measured thermistor resistance (ohms)
    !! output:
    !!         T:   temperature (C)
    !!

    !--------Arguments--------!
    REAL, INTENT(IN) :: rin
    REAL :: T

    !----------Local vars----------!
    REAL :: rth_log

    !--------Parameters------------!
    REAL(r8), PARAMETER :: a = 1.286212e-3
    REAL(r8), PARAMETER :: b = 2.355213e-4
    REAL(r8), PARAMETER :: c = 9.826046e-8
    REAL(r8), PARAMETER :: d = 8.835732e-8

    REAL, PARAMETER :: abs_zero = -273.16
    REAL, PARAMETER :: rmax = 4990.0

    IF (.NOT. Finite (rin)) THEN  ! Nothing to be done
       T = Qnan()
    ELSE IF (rin >= rmax) THEN  ! over the maximum ohms
       T = -100.0               ! minimum possible temperature
    ELSE IF (rin <= 0.0) THEN   ! under the minimum ohms
       T = 500.0                ! maximum possible temperature
    ELSE
       rth_log = LOG ((rmax * rin) / (rmax - rin)) ! log(Eq 5.5, pg 49)
       T = 1.0 / (a + rth_log * (b + rth_log * (c + rth_log * d))) + abs_zero
    ENDIF

  END FUNCTION Therm_temp

  FUNCTION PRD_temp (rin, r0) RESULT (T)

    !! See Eq. 5.4, Sec 5.1.1, pg 48 of MLS Level 1 ATB (v2.0)
    !!
    !!

    !! Convert PRD (ohms) to temperature (C)
    !! inputs:
    !!         rin: measured PRD resistance (ohms)
    !!         r0:  PRD resistance at 0 C (ohms)
    !! output:
    !!         T:   temperature (C)

    !--------Arguments--------!
    REAL, INTENT(IN) :: rin
    REAL, INTENT(IN) :: r0
    REAL :: T

    !----------Local vars----------!
    REAL(r8), PARAMETER :: a = 0.48945548411
    REAL(r8), PARAMETER :: b = 7.20107099888e-5
    REAL :: r

    r = rin * 500.0 / r0
    T = a * (r - 500.0) / (1.0 - b * r)

  END FUNCTION PRD_temp

  FUNCTION Use_equation (mnemonic, tval) RESULT (val)

    ! Use equation to scale analog telemetry

    ! inputs:
    !          mnemonic: telemetry mnemonic
    !          tval: input telemetry value

    !--------Arguments--------!
    CHARACTER (LEN=*), INTENT(IN) :: mnemonic
    REAL, INTENT(IN)  :: tval

    REAL :: val

    val = tval

    ! check mnemonic for equation to use

    IF (INDEX (mnemonic, "_Bus") /= 0) THEN
       IF (INDEX (mnemonic, "_V") /= 0) THEN
          IF (INDEX (mnemonic, "Quiet") /= 0) THEN
             val = (tval - 0.3113) / 0.1183   ! Quiet Bus Voltage
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = (tval - 0.2638) / 0.1192   ! Device Bus Voltage
          ENDIF
       ELSE IF (INDEX (mnemonic, "_I") /= 0) THEN
          IF (INDEX (mnemonic, "Quiet") /= 0) THEN
             val = (-0.156313 + SQRT (0.156313 * 0.156313 - 4.0 * 0.00026955 * &
                  (0.361124 - tval))) / (2 * 0.00026955)  ! Quiet Bus Current
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = 1.36 * (tval - 0.4852) / (1.3635 - 0.4852) ! Device Bus I
           ENDIF
       ENDIF
    ELSE IF (INDEX (mnemonic, "R4_IF_voltage") /= 0) THEN
       IF (tval > 0.0) THEN
          val = 8.4459 * LOG (tval) - 47.755  ! R4 IF Voltage
       ELSE
          val = QNan()
       ENDIF
    ENDIF

  END FUNCTION Use_equation


  SUBROUTINE ConvertEngCounts (GMAB_ON)

    !! Convert the raw counts into engineering units. The input engineering data
    !! is stored in the module variable EngPkt, which stores one MAF of
    !! Engineering data and is filled in each call to NextEngMaf, which calls
    !! this routine.

    USE EngTbls, ONLY: Eng_tbl, Riu_tbl, EngPkt, last_tlm_pt, temp_cal, &
     & vin_cal1, vin_cal2, prt1_cal1, prt1_cal2, prt2_cal1, prt2_cal2, &
     & ysi_cal1, ysi_cal2, Cal_const

    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info,MLSMSG_Warning
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC
    LOGICAL, DIMENSION(:) :: GMAB_ON

    INTEGER :: i, n, dn, pkt_no, riu_no, byte1,rs
    REAL :: scale

    REAL, PARAMETER :: abs_zero = -273.16
    REAL(r8) TAI93


    LOGICAL :: FoundBadPoints=.FALSE.

    CHARACTER(len=27) :: asciiUTC
    CHARACTER(len=256) :: msg
    INTEGER, EXTERNAL :: PGS_TD_EOSPMGIRDtoTAI
    logical, parameter :: DeeBug = .false.



! Extract the raw counts. For each RIU

    DO i = 1, SIZE(Riu_tbl) ! count RIUs
      
      ! Here 'i' is the RIU number.
      
      pkt_no = Riu_tbl(i)%pkt_no ! pkt for this RIU
      byte1 = Riu_tbl(i)%start_byte ! starting byte for this RIU

      Riu_tbl(i)%id_word = BigEndianStr (EngPkt(pkt_no)(byte1:byte1+1))

      DO n = Riu_tbl(i)%first_pt, Riu_tbl(i)%last_pt
        ! <whd> 
        ! {first,last}_pt are actually the first and last line numbers in the
        ! eng tables file (currently named engtlm.tbl) which mentions RIU(i);
        ! it is *not* a telemetry value!  Their purpose here is just to count
        ! the number of bytes that need to be extracted from this packet for
        ! this RIU

        ! So all that's really important here is the number of lines for this
        ! RIU, i.e.  nLines = Riu_tbl(i)%last_pt - Riu_tbl(i)%first_pt+1

        byte1 = byte1 + 2    ! next word

        ! so the packet is arranged as: (all quantities are 2-byte long)
        ! id 
        ! counts(mnem(1)) 
        ! counts(mnem(2))
        ! counts(mnem(3))
        ! ...
        ! counts(mnem(last_mnem)
        !!
        Eng_tbl(n)%counts = BigEndianStr (EngPkt(pkt_no)(byte1:byte1+1))


        ! <whd> why a special test for GM01_PRT1_Cal1 </whd>

        if ( n == 3 .and. Eng_tbl(n)%counts < 1 .and. DeeBug ) &
             & print *, 'Cal counts is <= 0 for Eng table at n ', n

        ! Save the calibration counts

        ! <whd> 
        ! cal_indx is index into the array EngTbls::Cal_Type_str, so
        ! Riu_tbl%cal_cnts has the same order as that array. Parameters that
        ! hold these indices are also defined in EngTbls ~line 49
        ! </whd>

        IF (Eng_tbl(n)%cal_indx > 0) THEN
          Riu_tbl(i)%Cal_cnts(Eng_tbl(n)%cal_indx) = Eng_tbl(n)%counts
          ! Don't know why ptr1_cal{1,2} also saves ptr2_cal{1,2}.
          IF (Eng_tbl(n)%cal_indx == prt1_cal1) THEN  ! save also as PRT-2's
            Riu_tbl(i)%Cal_cnts(prt2_cal1) = Eng_tbl(n)%counts
          ELSE IF (Eng_tbl(n)%cal_indx == prt1_cal2) THEN
            Riu_tbl(i)%Cal_cnts(prt2_cal2) = Eng_tbl(n)%counts
          ENDIF
        ENDIF

      ENDDO

    ENDDO

! Convert the counts

    ! <whd>
    ! `i' here is basically the lines of the engtlm.tbl file (excluding those in
    ! the `survival table', the file header, and section demarcating lines
    ! (i.e. those consisting of a line of dashes)
    ! </whd>

    DO i = 1, last_tlm_pt

      riu_no = Eng_tbl(i)%riu_no
      scale = Eng_tbl(i)%scale

      Eng_tbl(i)%value = QNan() ! assume the worst

      ! Skip the always bad THz sensor

      IF ((INDEX (Eng_tbl(i)%mnemonic, "THzAmbCalTgt_T2") /= 0) .OR. &
           (INDEX (Eng_tbl(i)%mnemonic, "THzAmbCalTgtRT2") /= 0)) CYCLE

      IF (riu_no >= 1 .AND. riu_no <= 4) THEN
        IF (.NOT. GMAB_ON(riu_no)) CYCLE     ! GMriu_no is NOT ON
        IF (MOD (riu_no, 2) == 0) THEN
          IF (GMAB_ON(riu_no-1)) CYCLE      ! Other related GMriu_no is ON
        ELSE
          IF (GMAB_ON(riu_no+1)) CYCLE      ! Other related GMriu_no is ON
        ENDIF
      ENDIF

      ! <whd>
      ! if id_word == 0, this RIU is off!
      ! Riu_tbl(riu_no)%id_word is the ID for the Riu_tbl(riu_no). The C&DH
      ! handbook (e.g. see pages 30-33, similarly for any fields marked
      ! [GMST]0[1-4]) Eng_tbl(i) the entry on the i-th non-trivial line in the
      ! engineering telemetry table file, currently called
      ! engtlm.tbl. Eng_tbl(i)%counts is read from the telemetry, using the pkt
      ! no and offset stored in Riu_tbl. Riu_tbl is filled in EngTbls.f90
      ! </whd>

      IF (Riu_tbl(riu_no)%id_word /= 0 .AND. &
           Eng_tbl(i)%counts > 0) THEN   
        ! RIU is on and there are non-zero counts: Data available to convert

        SELECT CASE (Eng_tbl(i)%type)

        CASE ("Cal")

          IF (Eng_tbl(i)%cal_indx /= temp_cal) THEN   ! not a temperature
            Eng_tbl(i)%value = Eng_tbl(i)%counts
          ELSE

            Eng_tbl(i)%value = Cal_freq ( Eng_tbl(i)%counts, & !measurement
                 & Riu_tbl(riu_no)%Cal_cnts(vin_cal1), & !high cal cnt
                 & Riu_tbl(riu_no)%Cal_cnts(vin_cal2), & !low cal cnt
                 & Cal_const(riu_no)%volts, & !high value
                 & 0.0) &                     !low value
                 & * 100.0 + abs_zero
          ENDIF

        CASE ("VIN")

          Eng_tbl(i)%value = Cal_freq (Eng_tbl(i)%counts, &
               & Riu_tbl(riu_no)%Cal_cnts(vin_cal1), &
               & Riu_tbl(riu_no)%Cal_cnts(vin_cal2), &
               & Cal_const(riu_no)%volts, 0.0)

          IF (scale /= 0.0) THEN
            Eng_tbl(i)%value = scale * Eng_tbl(i)%value
          ELSE
            Eng_tbl(i)%value = Use_equation (Eng_tbl(i)%mnemonic, &
                 & Eng_tbl(i)%value)
          ENDIF

        CASE ("PRT-1")

          dn = Eng_tbl(i)%counts
          IF (dn > 63000) dn = dn - 65536  ! adjust if DN is too high

          ! calls cal_freq inside PRT_temp. Commented to make that a bit
          ! clearer. Cases PRT-2, YSI-{1,2} are similar.

          Eng_tbl(i)%value = PRD_temp ( &
               & Cal_freq(dn, &
               &          Riu_tbl(riu_no)%Cal_cnts(prt1_cal1), &
               &          Riu_tbl(riu_no)%Cal_cnts(prt1_cal2), &
               &          Cal_const(riu_no)%prt1_hi, &
               &          Cal_const(riu_no)%prt1_low &
               &          ), &
               & scale)

        CASE ("PRT-2")

          Eng_tbl(i)%value = PRD_temp ( &
               Cal_freq( &
               &                       Eng_tbl(i)%counts, &
               &                       Riu_tbl(riu_no)%Cal_cnts(prt2_cal1), &
               &                       Riu_tbl(riu_no)%Cal_cnts(prt2_cal2), &
               &                       Cal_const(riu_no)%prt2_hi, &
               &                       Cal_const(riu_no)%prt2_low &
               &                       ), &
               &                scale)

        CASE ("YSI")

          Eng_tbl(i)%value = Therm_temp (&
               &              Cal_freq(&
               &                 Eng_tbl(i)%counts, &
               &                 Riu_tbl(riu_no)%Cal_cnts(ysi_cal1), &
               &                 Riu_tbl(riu_no)%Cal_cnts(ysi_cal2), &
               &                 Cal_const(riu_no)%therm_hi, Cal_const(riu_no)%therm_low &
               &                 ) &
               &               )

        CASE DEFAULT

          if ( DeeBug ) PRINT *, "Unknown tlm type!"  !! use standard routine here

        END SELECT

      ELSEIF ( i < 307 .or. i > 321 ) THEN
        !<whd:comment> 
!
        ! Riu_tbl(riu_no)%id_word /= 0 .OR.Eng_tbl(i)%counts <= 0, but i < 307
        ! or > 321. Lines 307-320 in engtlm.tbl correspond to GM15 and
        ! SM01. Comment below says those are always NaNs
!
        !</whd:comment> 
        if ( DeeBug ) PRINT *, 'Data not available to convert'
        IF (Eng_tbl(i)%mnemonic .eq. 'Spare') THEN
          if ( DeeBug ) print *,"But it's a `Spare', so we don't really care!"
        ELSE
          if ( DeeBug ) print *,"Mnemonic = ",Eng_tbl(i)%mnemonic
        ENDIF
      ENDIF
      ! Apparently values 307-320 are all NaNs Who wooda guessed!  
!
      ! Line numbers 307-320 in the engtlm.tbl file (not counting lines which
      ! exclude header, comments and section demarcating lines of all dashes)
      ! are for the GigaHertz switch network, which is only turned on and then
      ! immediately turned off, so these are, in a sense, bogus NaNs and are
      ! not reported. Hence their exclusion from the following code block
      IF ( isNaN(Eng_tbl(i)%value) .and. &
           & ( i < 307 .or. i > 321 ) &
           & ) THEN
        FoundBadPoints=.TRUE.
        if ( DeeBug ) then
          print *, i, 'th engineering value is a NaN'
          print *, 'type  ', Eng_tbl(i)%type
          print *, 'riu number ', riu_no
          print *, 'mnemonic ', trim(Eng_tbl(i)%mnemonic)
          if ( riu_no > 0 .and. riu_no < 5 ) print *, 'GMriu on? ', GMAB_ON(riu_no)
          if ( riu_no > 1 .and. riu_no < 6 .and. MOD (riu_no, 2) == 0 ) &
               & print *, 'prev GMriu on? ', GMAB_ON(riu_no-1)
          if ( riu_no > -1 .and. riu_no < 4 .and. MOD (riu_no, 2) /= 0 ) &
               & print *, 'next GMriu on? ', GMAB_ON(riu_no+1)
          print *, 'id_word ', Riu_tbl(riu_no)%id_word
          print *, 'counts ', Eng_tbl(i)%counts
        endif
      ENDIF
      ! So, it's because GMriu_no is NOT ON for these rius? Tell me, where is
      ! any of this documented? Does anyone here know anything? 
!
      ! (<whd> I'm presuming that vince wrote that last comment. As for me
      ! I have to say 'no', apparently no one here knows anything!</whd>)

    ENDDO ! loop over engtlm.tbl lines 

    IF (FoundBadPoints) THEN
       rs = PGS_TD_EOSPMGIRDtoTAI (engpkt(1)(8:15), TAI93)
       asciiUTC="BadTime"
       IF (rs /= 0) THEN 
          print *,"ConvertEngCounts: Can't get engineering packet TAI93 time!"
          CALL MLSMessage(MLSMSG_Warning,ModuleName, &
               &   "ConvertEngCounts: Can't get engineering packet TAI93 time!")          
       ELSE
          rs = PGS_TD_TAItoUTC (TAI93, asciiUTC)
       ENDIF
       write(msg,*) 'ConvertEngCounts: bad data! MIF time: '//trim(asciiUTC)
       print *,trim(msg)
       CALL MLSMessage(MLSMSG_Info,ModuleName, trim(msg))

    ENDIF
  END SUBROUTINE ConvertEngCounts

!=============================================================================
  SUBROUTINE NextEngMAF (more_data)
!=============================================================================

    USE EngTbls, ONLY: EngPkt, EngMAF, Eng_tbl
    USE L0Utils, ONLY: ReadL0Eng

    !! Get the next MAF's engineering data

    LOGICAL, INTENT (OUT) :: more_data

    LOGICAL :: data_OK, GMAB_ON(4)
    INTEGER :: maskbit3, maskbit7
    DATA maskbit3 / z'8' /
    DATA maskbit7 / z'80' /

! Save last A/B side readbacks ("A" side default):

    CHARACTER(len=1), SAVE :: ASE_Side = "A", GSM_Side = "A"

    DO
       CALL ReadL0Eng (EngPkt, EngMAF%MAFno, EngMAF%TotalMAF, &
            EngMAF%MIFsPerMAF, EngMAF%secTAI, data_OK, more_data)
       IF (.NOT. more_data) RETURN   ! nothing more to do
       IF (data_OK) EXIT
    ENDDO

    ! <whd:comment> 0.25 ~= 1.5 MIFs (MIF has a nominal dur of 1/6
    ! secs). Dominick thinks that the time in the engineering packets that goes
    ! with a particular set of science packets has a time tag about 1.5 MIFs
    ! later than the science packets they go with. However, in
    ! CalibWeightsFlags::ProcessMAFdata, the engineering data that's 'paired'
    ! with science data is selected from the previous Eng MAF.  See ~line 185,
    ! at the mention of NINT() </whd:comment>

    EngMAF%secTAI = EngMAF%secTAI - 0.25   ! vp's orig comment: actual time to line up with SCI

    ! Determine GM01 to GM04 A/B ON/OFF states:

    GMAB_ON(1) =(IAND (ICHAR (EngPkt(6)(103:103)), maskbit7) == maskbit7)
    GMAB_ON(2) =(IAND (ICHAR (EngPkt(6)(103:103)), maskbit3) == maskbit3)
    GMAB_ON(3) =(IAND (ICHAR (EngPkt(6)(109:109)), maskbit3) == maskbit3)
    GMAB_ON(4) =(IAND (ICHAR (EngPkt(6)(110:110)), maskbit7) == maskbit7)

    ! Determine which side ASE (GM01/GM02) is on:

    IF (GMAB_ON(1) .AND. .NOT. GMAB_ON(2)) THEN
       EngMAF%ASE_Side = "A"
    ELSE IF (.NOT. GMAB_ON(1) .AND. GMAB_ON(2)) THEN
       EngMAF%ASE_Side = "B"
    ELSE
       EngMAF%ASE_Side = ASE_Side    ! "U" - use previous known position instead
    ENDIF
    ASE_Side = EngMAF%ASE_Side

    ! Determine which side GSM (GM03/GM04) is on:

    IF (GMAB_ON(3) .AND. .NOT. GMAB_ON(4)) THEN
       EngMAF%GSM_Side = "A"
    ELSE IF (.NOT. GMAB_ON(3) .AND. GMAB_ON(4)) THEN
       EngMAF%GSM_Side = "B"
    ELSE
       EngMAF%GSM_Side = GSM_Side    ! "U" - use previous known position instead
    ENDIF
    GSM_Side = EngMAF%GSM_Side

    ! Convert engineering counts

    CALL ConvertEngCounts (GMAB_ON)

    !! Save the required data for later use:

    EngMAF%Eng%value = Eng_tbl%value
    EngMAF%Eng%mnemonic = Eng_tbl%mnemonic

  END SUBROUTINE NextEngMAF

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE EngUtils
!=============================================================================

! $Log$
! Revision 2.16  2024/10/10 20:15:54  pwagner
! Reduce routine messaging
!
! Revision 2.15  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.14.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.14  2015/01/21 19:30:19  pwagner
! Gets isNaN from MLSFillValues
!
! Revision 2.13  2015/01/14 00:32:53  pwagner
! Added debugging info for Calibration counts dropping to zero
!
! Revision 2.12  2008/05/06 15:56:44  perun
! Read eng packets until data good.
!
! Revision 2.11  2006/08/02 18:53:41  perun
! Define default A/B side reading and use if none available during reading
!
! Revision 2.10  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.9  2004/11/10 15:33:40  perun
! Adjust TAI time by -0.25 secs; change call to ReadL0Eng with additional flag
!
! Revision 2.8  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.7  2004/01/09 17:46:22  perun
! Version 1.4 commit
! Revision 2.6  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2001/02/23 18:55:17  perun
! Version 0.5 commit
!
