! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Options              !  Options and Settings for the MLSL2 program
!=============================================================================

  USE MLSMessageModule, only: MLSMSG_Error

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

  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  ! Set the following to TRUE before delivering level 2 to sips
  logical, private, parameter :: SIPS_VERSION =  .false. 
  logical            :: PUNISH_FOR_INVALID_PCF = SIPS_VERSION 
  logical, parameter :: PUNISH_FOR_NO_L1BRAD =   SIPS_VERSION
  logical, parameter :: PUNISH_FOR_NO_L1BOA =    SIPS_VERSION
  logical            :: PCF_FOR_INPUT =          SIPS_VERSION ! Open L2CF using PCF ?
  logical            :: PCF =                    SIPS_VERSION ! Use PCF ?
  logical            :: CREATEMETADATA =         SIPS_VERSION ! Create .met files ?
  logical            :: TOOLKIT =                SIPS_VERSION ! Use PGS_... routines?
  ! PCF controls whether the input and output file names are obtained
  ! from the PCF or the l2cf; if .false., the l2cf must supply every
  ! file name (L1B..) plus all the global settings (start, end times, ..)
  ! Note the following cascade of automatic negations:
  ! TOOLKIT=.false. =>  PCF=.false.
  ! PCF=.false.     =>  PCF_FOR_INPUT=.false.
  ! PCF=.false.     =>  PUNISH_FOR_INVALID_PCF=.false.
  ! PCF=.false.     =>  CREATEMETADATA=.false.
  ! PCF=.false.     =>  PENALTY_FOR_NO_METADATA=0

  ! Must files named in PCF have same case as short names used in l2cf?
  logical, parameter ::                          PCFL2CFSAMECASE = SIPS_VERSION
  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   ! Set the following to 1 before delivering to sips;
   ! when set to 0, it allows program to run w/o creating metadata
   integer            ::                         PENALTY_FOR_NO_METADATA = 0

   ! Set the following to -2 before delivering to sips;
   ! (its possible values and their effects on normal output:
   ! -1          sent to stdout (via print *, '...')
   ! -2          sent to Log file (via MLSMessage)
   ! < -2        both stdout and Log file
   ! > -1        Fortran 'unit=OUTPUT_PRINT_UNIT')
   integer            :: OUTPUT_PRINT_UNIT = -2

   ! Set the following to MLSMSG_Error before delivering to sips;
   ! when set higher, it allows program keep going despite errors
   ! when set lower, the program would quit even on warnings
   integer, parameter :: QUIT_ERROR_THRESHOLD = MLSMSG_Error

   ! Set the following to 2 before delivering to sips;
   ! If 0, program simply stops both upon normal termination
   ! as well as some abnormal ones (e.g. in parser)
   ! if 2, status will be 2 only if run complete
   ! and without error
   integer, parameter :: NORMAL_EXIT_STATUS = 2

!=============================================================================
END MODULE MLSL2Options
!=============================================================================

!
! $Log$
! Revision 2.11  2001/09/28 17:57:47  pwagner
! SIPS_VERSION controls other logical options
!
! Revision 2.10  2001/07/16 23:43:15  pwagner
! With settable NORMAL_EXIT_STATUS
!
! Revision 2.9  2001/05/30 22:56:48  pwagner
! Moved PCFL2CFSAMECASE here from OutputAndClose
!
! Revision 2.8  2001/05/15 23:46:08  pwagner
! Removed 2 settings from MLSL2Opts; now in switches
!
! Revision 2.7  2001/05/11 23:48:23  pwagner
! Changed to not echo globals; added note on SIPS
!
! Revision 2.6  2001/05/09 23:34:13  pwagner
! Added ECHO_GLOBAL_STNGS LOG_TO_STDOUT
!
! Revision 2.5  2001/05/06 20:54:40  pwagner
! Default settings should work for most jpl users
!
! Revision 2.4  2001/05/04 22:54:31  pwagner
! Added TOOLKIT, CREATEMETADATA, PCF_FOR_INPUT
!
! Revision 2.3  2001/04/20 20:41:52  pwagner
! Added QUIT_ERROR_THRESHOLD
!
! Revision 2.2  2001/04/17 20:26:28  pwagner
! Added OUTPUT_PRINT_UNIT
!
! Revision 2.1  2001/04/16 23:53:10  pwagner
! First commit
!
