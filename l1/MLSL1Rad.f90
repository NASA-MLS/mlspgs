! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Rad              ! Radiance data types for the MLSL1 program
!=============================================================================

  USE MLSL1Common
  USE MLSSignalNomenclature
  USE MLSMessageModule
  USE MLSStrings

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
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

  TYPE (Radiance_T), DIMENSION(:), POINTER :: L1Brad, FBrad, MBrad, WFrad

CONTAINS

!=============================================================================
  SUBROUTINE InitRad
!=============================================================================

    USE MLSL1Config, ONLY: MIFsGHz, MIFsTHz

    ! Arguments

    ! Local

    INTEGER :: i, status,RADMIFs
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal => NULL()
    CHARACTER (LEN=3) :: bandNo

    ! For this version, allocate for only standard and mid-band filter banks

    ALLOCATE (L1Brad(FBnum+MBnum+WFnum), STAT=status)
    IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
         & MLSMSG_Allocate//"L1BRad")

    ! Allocate and initialize for the standard filter banks

    DO i = 1, FBnum

       IF (i < 15) THEN
          RADMIFs = MIFsGHz
       ELSE
          RADMIFs = MIFsTHz
       ENDIF
       ALLOCATE (L1Brad(i)%value(FBchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBvalue")
       ALLOCATE (L1Brad(i)%precision(FBchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBprecision")
       WRITE (bandNo, '("B",i2.2)') i
       CALL ParseMLSSignalRequest (bandNo, signal, .false.)
       L1Brad(i)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    RADMIFs = MIFsGHz          ! Size for remainder of channels

    FBrad => L1Brad(1:FBnum)   ! Point to FB data

    ! Allocate and initialize for the mid-band filter banks

    DO i = 1, MBnum

       ALLOCATE (L1BRad(i+FBnum)%value(MBchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBvalue")
       ALLOCATE (L1BRad(i+FBnum)%precision(MBchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBprecision")
       WRITE (bandNo, '("B",i2.2)') (i+26)   ! start at band 27
       CALL ParseMLSSignalRequest (bandNo, signal, .false.)
       L1Brad(i+FBnum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    MBrad => L1Brad(FBnum+1:FBnum+MBnum)  ! Point to MB data

    ! Allocate and initialize for the wide filter banks

    DO i = 1, WFnum

       ALLOCATE (L1BRad(i+FBnum+MBnum)%value(WFchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"WFvalue")
       ALLOCATE (L1BRad(i+FBnum+MBnum)%precision(WFchans,RADMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"WFprecision")
       WRITE (bandNo, '("B",i2.2)') (i+31)   ! start at band 32
       CALL ParseMLSSignalRequest (bandNo, signal, .false.)
       L1Brad(i+FBnum+MBnum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    WFrad => L1Brad(FBnum+MBnum+1:FBnum+MBnum+WFnum)  ! Point to WF data

    RETURN

  END Subroutine InitRad

!=============================================================================
END MODULE MLSL1Rad
!=============================================================================

!
! $Log$
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/16 21:41:13  perun
! Added calls to MLSSignals routine
!

