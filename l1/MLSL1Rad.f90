! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
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

CONTAINS

!=============================================================================
  SUBROUTINE InitRad (L1Brad)
!=============================================================================

    ! Arguments

    TYPE (Radiance_T), DIMENSION(:), POINTER :: L1Brad

    ! Local

    INTEGER :: i, status
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal
    CHARACTER (LEN=3) :: bandNo

    ! For this version, allocate for only standard and mid-band filter banks

    ALLOCATE (L1Brad(FBnum+MBnum), STAT=status)
    IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
         & MLSMSG_Allocate//"L1BRad")

    ! Allocate and initialize for the standard filter banks

    DO i = 1, FBnum

       ALLOCATE (L1Brad(i)%value(FBchans,MaxMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBvalue")
       ALLOCATE (L1Brad(i)%precision(FBchans,MaxMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"FBprecision")
       WRITE (bandNo, '("B",i2.2)') i
       CALL ParseMLSSignalRequest (bandNo, signal, .false.)
       L1Brad(i)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    ! Allocate and initialize for the mid-band filter banks

    DO i = 1, MBnum

       ALLOCATE (L1BRad(i+FBnum)%value(MBchans,MaxMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBvalue")
       ALLOCATE (L1BRad(i+FBnum)%precision(MBchans,MaxMIFs), STAT=status)
       IF (status /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName,&
            & MLSMSG_Allocate//"MBprecision")
       WRITE (bandNo, '("B",i2.2)') (i+26)   ! start at band 27
       CALL ParseMLSSignalRequest (bandNo, signal, .false.)
       L1Brad(i+FBnum)%signal = signal(1)
       DEALLOCATE (signal)

    END DO

    RETURN

  END Subroutine InitRad

!=============================================================================
END MODULE MLSL1Rad
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.2  2000/02/16 21:41:13  perun
! Added calls to MLSSignals routine
!

