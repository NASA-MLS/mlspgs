! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE IEEE_ARITHMETIC              ! Common utilities for the MLSL1 program
!=============================================================================

!  use SUNSOWNIEEE, only: ir_isnan, ir_isinf, r_quiet_nan

  IMPLICIT NONE

  ! These are the ieee functions and constants needed in mls
  ! that Sun fails to supply.
  ! An improvement might be to introduce function interfaces to allow
  ! both single and double precision versions.
  public :: ieee_value, ieee_quiet_nan, ieee_is_finite

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSFile:$"
  !---------------------------------------------------------------------------

  ! This module contains glue routines for Sun's own f95 compiler
  ! since it fails to conform to ISO/IEC TR15580:1998(E).
  ! If we should ever obtain one that conforms we may cheerfully
  ! dispose of this crude hack

  private

  ! The following is made necessary only because Sun's own f95
  ! fails to comply--this obviously fails if any number ever matches this

  real, parameter :: ieee_quiet_nan = -9.E19
!  real, parameter :: ieee_quiet_nan = r_quiet_nan()

CONTAINS

  FUNCTION IEEE_VALUE( ARG1, ARG2 ) result ( the_value )
  ! Formal args
    real, intent(in) ::          arg1
    real, intent(in) ::          arg2
    real ::                      the_value
    
  ! Private
    the_value = arg2
  END FUNCTION IEEE_VALUE
  
  LOGICAL FUNCTION IEEE_IS_FINITE( ARG )
  ! Formal args
    real, intent(in) ::          arg
!    integer, external ::         ir_isnan
!    integer, external ::         ir_isinf
  ! Private
    
    IEEE_IS_FINITE = .FALSE.
!    if( ir_isnan(arg) == 1 ) RETURN
!    if( ir_isinf(arg) == 1 ) RETURN
    IEEE_IS_FINITE = .TRUE.
  END FUNCTION IEEE_IS_FINITE
  
END MODULE IEEE_ARITHMETIC

!
! $Log$
! Revision 1.2  2001/09/14 17:35:33  pwagner
! Fixed an obvious blunder; expanded some comments
!
! Revision 1.1  2001/09/13 23:18:37  pwagner
! First commit
!
