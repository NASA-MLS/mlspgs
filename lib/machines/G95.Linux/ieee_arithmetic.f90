! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

!=============================================================================
MODULE IEEE_ARITHMETIC              ! Common utilities for the MLSL1 program
!=============================================================================

!  use g95SOWNIEEE, only: ir_isnan, ir_isinf, r_quiet_nan

  implicit NONE
  private

  ! These are the ieee functions and constants needed in MLS
  ! that g95 fails to supply.

  public :: IEEE_Class_Type
  public :: IEEE_Is_finite, IEEE_Is_NaN, IEEE_Quiet_NaN, IEEE_Support_NaN
  public :: IEEE_Value

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module contains glue routines for g95 compiler
  ! since it fails to conform to ISO/IEC TR15580:1998(E).
  ! If we should ever obtain one that conforms we may cheerfully
  ! dispose of this crude hack

  type IEEE_Class_Type
    private
    integer :: What        ! 1 = IEEE_Quiet_Nam
                           ! 2 = IEEE_Signaling_Nan
  end type IEEE_Class_Type

  type(ieee_class_type), parameter :: ieee_quiet_nan = ieee_class_type(1)
  type(ieee_class_type), parameter :: ieee_signaling_nan = ieee_class_type(2)

!  real, parameter :: ieee_quiet_nan = r_quiet_nan()

  interface IEEE_Is_NaN
    module procedure IEEE_Is_NaN_D, IEEE_Is_NaN_S
  end interface

  interface IEEE_Support_NaN
    module procedure IEEE_Support_NaN_All
    module procedure IEEE_Support_NaN_D, IEEE_Support_NaN_S
  end interface

  interface IEEE_IS_FINITE
    module procedure IEEE_IS_FINITE_D, IEEE_IS_FINITE_S
  end interface

  interface IEEE_Value
    module procedure IEEE_Value_D, IEEE_Value_S
  end interface

CONTAINS

  LOGICAL FUNCTION IEEE_IS_FINITE_S( ARG )
  ! Formal args
    real, intent(in) ::          arg
!    integer, external ::         ir_isnan
!    integer, external ::         ir_isinf
  ! Private
    
    IEEE_IS_FINITE_S = .FALSE.
    if( IEEE_Is_NaN_S(arg) ) RETURN
    if( IEEE_Is_Inf_io_S(arg) ) RETURN
!    if( ir_isnan(arg) == 1 ) RETURN
!    if( ir_isinf(arg) == 1 ) RETURN
    IEEE_IS_FINITE_S = .TRUE.
  END FUNCTION IEEE_IS_FINITE_S
  
  LOGICAL FUNCTION IEEE_IS_FINITE_D( ARG )
  ! Formal args
    double precision, intent(in) ::          arg
!    integer, external ::         ir_isnan
!    integer, external ::         ir_isinf
  ! Private
    
    IEEE_IS_FINITE_D = .FALSE.
    if( IEEE_Is_NaN_D(arg) ) RETURN
    if( IEEE_Is_Inf_io_D(arg) ) RETURN
!    if( ir_isnan(arg) == 1 ) RETURN
!    if( ir_isinf(arg) == 1 ) RETURN
    IEEE_IS_FINITE_D = .TRUE.
  END FUNCTION IEEE_IS_FINITE_D
  
  logical function IEEE_Is_Inf_io_D ( X ) result(res)
    double precision, intent(in) :: X
    character(len=80) :: reschar
    write(reschar, *) x
    res = ( index(reschar, 'Inf') > 0 .or. &
      & index(reschar, 'inf') > 0 .or. index(reschar, 'INF') > 0 )
  end function IEEE_Is_Inf_io_D

  logical function IEEE_Is_Inf_io_S ( X ) result(res)
    real, intent(in) :: X
    character(len=80) :: reschar
    write(reschar, *) x
    res = ( index(reschar, 'Inf') > 0 .or. &
      & index(reschar, 'inf') > 0 .or. index(reschar, 'INF') > 0 )
  end function IEEE_Is_Inf_io_S

  elemental logical function IEEE_Is_NaN_D ( X )
    double precision, intent(in) :: X
    IEEE_Is_NaN_D = .not. ( x <= 0.0 .or. x >= 0.0 )
  end function IEEE_Is_NaN_D

  elemental logical function IEEE_Is_NaN_S ( X )
    real, intent(in) :: X
    IEEE_Is_NaN_S = .not. ( x <= 0.0 .or. x >= 0.0 )
  end function IEEE_Is_NaN_S

  logical function IEEE_Support_NaN_All ( )
    IEEE_Support_NaN_All = .true.
  end function IEEE_Support_NaN_All

  logical function IEEE_Support_NaN_D ( X )
    double precision, intent(in) :: X
    IEEE_Support_NaN_D = .true.
  end function IEEE_Support_NaN_D

  logical function IEEE_Support_NaN_S ( X )
    real, intent(in) :: X
    IEEE_Support_NaN_S = .true.
  end function IEEE_Support_NaN_S

  FUNCTION IEEE_VALUE_D ( X, CLASS ) result ( the_value )
  ! Formal args
    double precision, intent(in) ::      X
    type(ieee_class_type), intent(in) :: CLASS
    double precision ::                  the_value
    
  ! Private

  ! The following is made necessary only because g95
  ! fails to comply--this obviously fails if any number ever matches this
    select case ( class%what )
    case ( ieee_quiet_nan%what )     ! IEEE_Quiet_NaN
      the_value = -9.e19
    case ( ieee_signaling_nan%what ) ! IEEE_Signaling_NaN
      the_value = 9.e19
    end select
  END FUNCTION IEEE_VALUE_D

  FUNCTION IEEE_VALUE_S ( X, CLASS ) result ( the_value )
  ! Formal args
    real, intent(in) ::                  X
    type(ieee_class_type), intent(in) :: CLASS
    real ::                              the_value
    
  ! Private

  ! The following is made necessary only because g95
  ! fails to comply--this obviously fails if any number ever matches this
    select case ( class%what )
    case ( ieee_quiet_nan%what )     ! IEEE_Quiet_NaN
      the_value = -9.e19
    case ( ieee_signaling_nan%what ) ! IEEE_Signaling_NaN
      the_value = 9.e19
    end select
  END FUNCTION IEEE_VALUE_S

END MODULE IEEE_ARITHMETIC

!
! $Log$
