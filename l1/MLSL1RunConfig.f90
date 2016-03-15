! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSL1RunConfig! Stores run configuration
!=============================================================================

  ! This module is meant to store whatever information about the configuration
  ! of the run might be needed down stream. At the moment (Fri Aug 7 2015) the
  ! only thing I need is the name of the executable, which must be set in each
  ! of the `mains' (mlsl0sn, mlsl1log, mlsl1g and mlsl1t), and it primarily for
  ! use with the MLSL1Debug and SnoopMLSL1 modules.
  !
  ! But who knows? Maybe it will become more useful later.

  IMPLICIT NONE
  PRIVATE
  
  
  CHARACTER (len=256) :: MLSL1Executable
  PUBLIC MLSL1Executable

END MODULE MLSL1RunConfig
! $Id$
!
! Modifications
!
! $Log$
! Revision 1.1.2.1  2016/02/25 19:50:37  whdaffer
! Initial revision
!
!
