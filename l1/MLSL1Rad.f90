! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Rad     ! Radiance data types and routines for the MLSL1 program
!=============================================================================

  USE MLSL1Common, ONLY: FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, R4, &
       DACSnum, DACSchans, GHzNum, THzNum, THzChans, Chan_R_T
  USE MLSSignalNomenclature !, ONLY: ParseMLSSignalRequest, MLSSignal_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Allocate

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: InitRad, BandToBanks, UpdateRadSignals, RadPwr
  PUBLIC :: Radiance_T, L1Brad, FBrad, MBrad, WFrad, DACSrad, THzRad

  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module defines the radiance data types and the nominal
  ! nomenclature sets for the MLSL1 program.

  !---------------------------------------------------------------------------

  TYPE Radiance_T
     TYPE (MLSSignal_T) :: signal
     REAL(r4), DIMENSION (:,:), POINTER :: value, precision
  END TYPE Radiance_T

  TYPE (Radiance_T), DIMENSION(:), POINTER :: L1Brad, FBrad, MBrad, WFrad, &
       DACSrad, THzRad

CONTAINS

!=============================================================================
  SUBROUTINE InitRad (THz)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, MIFsGHz, MIFsTHz
    USE MLSL1Common, ONLY: BandSwitch

    ! Arguments

    LOGICAL :: THz

    ! Local

    INTEGER :: i, status, BandNo, BankNo
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal => NULL()
    CHARACTER (LEN=11) :: request
    INTEGER, PARAMETER :: DACSbandNo(4) = (/ 25, 23, 24, 22 /)

    IF (THz) THEN  ! Allocate for THz
       ALLOCATE (L1Brad(THzNum), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"L1BRad")
       THzRad => L1Brad
       DO i = 1, THzNum

          ALLOCATE (L1Brad(i)%value(THzchans,MIFsTHz), STAT=status)
          IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
               & MLSMSG_Allocate//"FBvalue")
          ALLOCATE (L1Brad(i)%precision(THzchans,MIFsTHz), STAT=status)
          IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
               & MLSMSG_Allocate//"FBprecision")

          BankNo = i + 14   !start at 15
          SELECT CASE (i)
          CASE (1)
             BandNo = BandSwitch(5)
          CASE (THzNum)
             BandNo = BandSwitch(4)
             BankNo = 12
          CASE DEFAULT
             BandNo = BankNo
          END SELECT

          WRITE (request, '("B",i2.2,".FB25-",i2)') BandNo, BankNo
          request = request(1:9)//ADJUSTL(request(10:11))
          CALL ParseMLSSignalRequest (request, signal, .FALSE.)

          L1Brad(i)%signal = signal(1)
          DEALLOCATE (signal)

       END DO

       RETURN  ! Nothing more to do

    ENDIF

    ALLOCATE (L1Brad(GHzNum+MBnum+WFnum+DACSnum), STAT=status)
    IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
         & MLSMSG_Allocate//"L1BRad")

    ! Allocate and initialize for the standard filter banks

    DO i = 1, GHzNum

       ALLOCATE (L1Brad(i)%value(FBchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBvalue")
       ALLOCATE (L1Brad(i)%precision(FBchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBprecision")

       SELECT CASE (i)
       CASE (3)
          BandNo = BandSwitch(2)
       CASE (8)
          BandNo = BandSwitch(3)
       CASE (12)
          BandNo = BandSwitch(4)
       CASE DEFAULT
          BandNo = i
       END SELECT

       WRITE (request, '("B",i2.2,".FB25-",i2)') BandNo, i
       request = request(1:9)//ADJUSTL(request(10:11))
       CALL ParseMLSSignalRequest (request, signal, .FALSE.)
       L1Brad(i)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    FBrad => L1Brad(1:GHzNum)   ! Point to FB data

    ! Allocate and initialize for the mid-band filter banks

    DO i = 1, MBnum

       ALLOCATE (L1BRad(i+GHzNum)%value(MBchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBvalue")
       ALLOCATE (L1BRad(i+GHzNum)%precision(MBchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBprecision")
       WRITE (request, '("B",i2.2,".MB11-",i1)') (i+26), i
       CALL ParseMLSSignalRequest (request, signal, .FALSE.)
       L1Brad(i+GHzNum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    MBrad => L1Brad(GHzNum+1:GHzNum+MBnum)  ! Point to MB data

    ! Allocate and initialize for the wide filter banks

    DO i = 1, WFnum

       ALLOCATE (L1BRad(i+GHzNum+MBnum)%value(WFchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"WFvalue")
       ALLOCATE (L1BRad(i+GHzNum+MBnum)%precision(WFchans,MIFsGHz), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"WFprecision")
       WRITE (request, '("B",i2.2,".WF4-",i1)') (i+31), i
       CALL ParseMLSSignalRequest (request, signal, .FALSE.)
       L1Brad(i+GHzNum+MBnum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    WFrad => L1Brad(GHzNum+MBnum+1:GHzNum+MBnum+WFnum)  ! Point to WF data

    ! Allocate and initialize for the DACS filter banks

    DO i = 1, DACSnum

       ALLOCATE (L1BRad(i+GHzNum+MBnum+WFnum)%value(DACSchans,MIFsGHz), &
            STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"DACSvalue")
       ALLOCATE (L1BRad(i+GHzNum+MBnum+WFnum)%precision(DACSchans,MIFsGHz), &
            STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"DACSprecision")
       IF (i == 1) THEN
          BandNo = BandSwitch(1)
       ELSE
          BandNo = DACSbandNo(i)
       ENDIF
       WRITE (request, '("B",i2.2,".DACS-",i1)') BandNo, i
       CALL ParseMLSSignalRequest (request, signal, .FALSE.)
       L1Brad(i+GHzNum+MBnum+WFnum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    DACSrad => L1Brad(GHzNum+MBnum+WFnum+1:GHzNum+MBnum+WFnum+DACSnum)

  END SUBROUTINE InitRad

!=============================================================================
  SUBROUTINE UpdateRadSignals (BandSwitch)
!=============================================================================

    USE MLSL1Common, ONLY: SwitchBank

    INTEGER :: BandSwitch(*)

    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal => NULL()
    CHARACTER (LEN=11) :: request
    INTEGER :: i

!! THz first

    IF (SIZE(L1Brad) == THzNum) THEN
       DO i = 4, 5   ! Switches 4 & 5
          IF (BandSwitch(i) > 0) THEN
             WRITE (request, '("B",i2.2,".FB25-",i2)') BandSwitch(i), &
                  SwitchBank(i)
             request = request(1:9)//ADJUSTL(request(10:11))
             CALL ParseMLSSignalRequest (request, signal, .FALSE.)
             IF (i == 5) THEN
                L1Brad(1)%signal = signal(1)
             ELSE
                L1Brad(THzNum)%signal = signal(1)
             ENDIF
             DEALLOCATE (signal)
          ENDIF
       ENDDO

       RETURN    ! Nothing more

    ENDIF

!! DACS

    IF (BandSwitch(1) > 0) THEN
       WRITE (request, '("B",i2.2,".DACS-1")') BandSwitch(1)
       CALL ParseMLSSignalRequest (request, signal, .FALSE.)
       L1Brad(1+GHzNum+MBnum+WFnum)%signal = signal(1)
       DEALLOCATE (signal)
    ENDIF
    
!! FB's next

    DO i = 2, 4
       IF (BandSwitch(i) > 0) THEN
          WRITE (request, '("B",i2.2,".FB25-",i2)') BandSwitch(i), SwitchBank(i)
          request = request(1:9)//ADJUSTL(request(10:11))
          CALL ParseMLSSignalRequest (request, signal, .FALSE.)
          L1Brad(SwitchBank(i))%signal = signal(1)
          DEALLOCATE (signal)
       ENDIF
    ENDDO

  END SUBROUTINE UpdateRadSignals

  SUBROUTINE BandToBanks (band, bank)

    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal => NULL()
    INTEGER :: band, bank(2)
    CHARACTER(LEN=3) :: request

    bank = 0   ! nothing yet
    WRITE (request, '("B", i2.2)') band
    CALL ParseMLSSignalRequest (request, signal, .FALSE.)
    bank(1) = minval (signal%spectrometernumber)
    bank(2) = maxval (signal%spectrometernumber)

    DEALLOCATE (signal)

  END SUBROUTINE BandToBanks

  FUNCTION RadPwr (Hz, T) RESULT (P)

    USE MLSL1Common, ONLY: boltz, planck

!! Calculate radiant power per unit bandwidth

    REAL :: Hz    !! frequency in Hz
    REAL :: T     !! temperature in Kelvin
    REAL :: P     !! radiant power

    P = (planck * Hz) / (boltz * (exp ((planck * Hz) / (boltz * T)) - 1.0))

  END FUNCTION RadPwr

!=============================================================================
END MODULE MLSL1Rad
!=============================================================================

!
! $Log$
! Revision 2.3  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/16 21:41:13  perun
! Added calls to MLSSignals routine
!

