! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Options              !  Options and Settings for the MLSL2 program
!=============================================================================

  USE MLSCommon
  USE MLSMessageModule

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module simply contains initial or permanent settings. Values
  ! are chosen according to what is most suitable for the environment.
  ! For example
  ! certain settings may be appropriate during development but not
  ! for production use, i.e. sips. Therefore, this is a convenient place
  ! to hold everything that needs to be changed before delivery.
  
  ! See also MLSL2Common and MLSL2PCF

  ! --------------------------------------------------------------------------

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Set the following to TRUE before delivering level 2 to sips
  logical, parameter :: PUNISH_FOR_INVALID_PCF=.false.  ! set to true
  logical, parameter :: PUNISH_FOR_NO_L1BRAD=.false.  ! set to true
  logical, parameter :: PUNISH_FOR_NO_L1BOA=.false.  ! set to true
  logical :: PCF = .false.                         ! Open L2CF using PCF ?
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   ! Change the following to 1 before delivering to sips;
   ! when set to 0, it allows program to run w/o creating metadata
   integer, parameter :: PENALTY_FOR_NO_METADATA = 0

   ! Change the following to -2 before delivering to sips;
   ! (its possible values and their effects on normal output:
   ! -1          sent to stdout (via print *, '...')
   ! -2          sent to Log file (via MLSMessage)
   ! < -2        both stdout and Log file
   ! > -1        Fortran 'unit=OUTPUT_PRINT_UNIT')
   integer, parameter :: OUTPUT_PRINT_UNIT = -1

!=============================================================================
END MODULE MLSL2Options
!=============================================================================

!
! $Log$
! Revision 2.2  2001/04/17 20:26:28  pwagner
! Added OUTPUT_PRINT_UNIT
!
! Revision 2.1  2001/04/16 23:53:10  pwagner
! First commit
!
