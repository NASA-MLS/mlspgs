! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Radiances ! Determine radiances
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
       MaxMIFs, absZero_C, BandWidth, tau
  USE MLSL1Utils, ONLY : GetIndexedAvg
  USE EngTbls, ONLY : Eng_MAF_T, CalTgtIndx
  USE Calibration, ONLY : CalWin, CalMAFdata, limb_time, limb_counts, &
       space_interp, target_interp, space_err, target_err, space_weight, &
       target_weight, zero_counts
  USE MLSL1Rad, ONLY : FBrad, MBrad, WFrad

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE CalcRadiance (limb_counts, space_counts, target_counts, &
       zero_counts, sum_w2, sum_wg2, space_temp, target_temp, bandwidth, &
       rad, rad_err)

    REAL(r8) :: limb_counts, space_counts, target_counts, zero_counts
    REAL(r8) :: sum_w2, sum_wg2
    REAL :: space_temp, target_temp, bandwidth, rad, rad_err

    REAL :: gain

    rad = 0.0        ! nothing yet
    rad_err = -1.0   ! nothing yet

    IF (limb_counts == 0.0d0 .OR. space_counts == target_counts) RETURN

    gain =  (target_counts - space_counts) / (target_temp - space_temp)
    rad =  (limb_counts - space_counts) / gain + space_temp

    rad_err = (limb_counts - zero_counts)**2 + (space_counts - zero_counts)**2 &
         * sum_w2 + (limb_counts - space_counts)**2 * &
         ((target_counts - zero_counts) / (target_counts - space_counts))**2 &
         * (1.0 + sum_w2) * sum_wg2
    rad_err = sqrt (rad_err / (bandwidth * tau)) / gain

  END SUBROUTINE CalcRadiance

  SUBROUTINE CalcLimbRads

    USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz

    TYPE (Eng_MAF_T) :: EMAF
    TYPE ValPrec_T
       REAL :: value
       REAL :: precision
    END TYPE ValPrec_T
    TYPE Rad_T
       TYPE (ValPrec_T) :: FB(FBchans,FBnum)    ! standard filter banks
       TYPE (ValPrec_T) :: MB(MBchans,MBnum)    ! mid-band filter banks
       TYPE (ValPrec_T) :: WF(WFchans,WFnum)    ! wide filters
    END TYPE Rad_T
    TYPE (Rad_T) :: Rad(0:MaxMIFs-1)

    INTEGER :: i, j, bank, chan, MIF_index, MIF_index_MAX
    INTEGER :: time_index, start_index, end_index, windex, CalWin_end
    REAL :: GHz_T1, GHz_T2, THz_T1, target_T
    REAL :: space_T
    REAL(r8) :: C_zero

    EMAF = CalWin%MAFdata(CalWin%central)%EMAF

    GHz_T1 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzAmb) + absZero_C
!!$    print *, "GHzAmb avg: ", GHz_T1
    GHz_T2 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%GHzCntl) + absZero_C
!!$    print *, "GHzCntl avg: ", GHz_T2
    THz_T1 = GetIndexedAvg (EMAF%eng%value, CalTgtIndx%THzAmb) + absZero_C
!!$    print *, "THzAmb avg: ", THz_T1

    GHz_T1 = 300.0    !!! For simulation TEST
    space_T = 0.0     !!! TEST FOR NOW!!!

    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    CalWin_end = CalWin%MAFdata(CalWin%size)%end_index

    DO time_index = start_index, end_index   ! for every MIF in the MAF

       MIF_index = time_index - start_index  ! MIF # within the central MAF

! Interpolate limb values

       DO bank = 1, FBnum

! Determine which Target temp to use

          target_T = GHz_T1

! Determine MAX allowable MIF index

          IF (bank < 15) THEN
             MIF_index_MAX = MIFsGHz - 1
          ELSE
             MIF_index_MAX = MIFsTHz - 1
          ENDIF

          IF (MIF_index <= MIF_index_MAX) THEN
             DO chan = 1, FBchans
          !C_zero = MAXVAL (zero_counts(start_index:end_index)%FB(chan,bank))
                C_zero = 1500
                CALL CalcRadiance (limb_counts(time_index)%FB(chan,bank), &
                     space_interp(MIF_index)%FB(chan,bank), &
                     target_interp(MIF_index)%FB(chan,bank), C_zero, &
                     space_err(MIF_index)%FB(chan,bank), &
                     target_err(MIF_index)%FB(chan,bank), &
                     space_T, target_T, BandWidth%FB(chan,bank), &
                     FBrad(bank)%value(chan,MIF_index+1), &
                     FBrad(bank)%precision(chan,MIF_index+1))
             ENDDO
          ENDIF
       ENDDO

       MIF_index_MAX = MIFsGHz - 1

       IF (MIF_index <= MIF_index_MAX) THEN
          DO bank = 1, MBnum
             DO chan = 1, MBchans
             !C_zero = MAXVAL (zero_counts(start_index:end_index)%MB(chan,bank))
                C_zero = 1500
                CALL CalcRadiance (limb_counts(time_index)%FB(chan,bank), &
                     space_interp(MIF_index)%MB(chan,bank), &
                     target_interp(MIF_index)%MB(chan,bank), C_zero, &
                     space_err(MIF_index)%MB(chan,bank), &
                     target_err(MIF_index)%MB(chan,bank), &
                     space_T, target_T, BandWidth%MB(chan,bank), &
                     MBrad(bank)%value(chan,MIF_index+1), &
                     MBrad(bank)%precision(chan,MIF_index+1))
             ENDDO
          ENDDO

          DO bank = 1, WFnum
             DO chan = 1, WFchans
             !C_zero = MAXVAL (zero_counts(start_index:end_index)%WF(chan,bank))
                C_zero = 1500
                CALL CalcRadiance (limb_counts(time_index)%WF(chan,bank), &
                     space_interp(MIF_index)%WF(chan,bank), &
                     target_interp(MIF_index)%WF(chan,bank), C_zero, &
                     space_err(MIF_index)%WF(chan,bank), &
                     target_err(MIF_index)%WF(chan,bank), &
                     space_T, target_T, BandWidth%WF(chan,bank), &
                     WFrad(bank)%value(chan,MIF_index+1), &
                     WFrad(bank)%precision(chan,MIF_index+1))
             ENDDO
          ENDDO
       ENDIF
    ENDDO

!!$DO i = 0, MIF_index_MAX
!!$  print *, INT(limb_counts(start_index+i)%FB(1,1)), &
!!$       INT(space_interp(i)%FB(1,1)),&
!!$       INT(target_interp(i)%FB(1,1)), FBrad(1)%value(1,i+1), &
!!$       FBrad(1)%precision(1,i+1)
!!$ENDDO

  END SUBROUTINE CalcLimbRads

!=============================================================================
END MODULE Radiances
!=============================================================================

! $Log$
! Revision 2.2  2001/09/10 16:17:56  perun
! Added CalMAFdata from Calibration module
!
! Revision 2.1  2001/02/23 20:55:04  perun
! Version 0.5 commit
!
