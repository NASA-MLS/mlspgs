! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE IEEE_ARITHMETIC              ! Common utilities for the MLSL1 program
!=============================================================================

  use Inf_NaN_Detection, only: isnan, isinf, isposinf, isneginf

  implicit NONE
  private

  ! These are the ieee functions and constants needed in MLS
  ! that Lahey fails to supply.

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

  ! This module contains glue routines for Lahey's own f95 compiler
  ! since it fails to conform to ISO/IEC TR15580:1998(E).
  ! If we should ever obtain one that conforms we may cheerfully
  ! dispose of this crude hack


  type IEEE_Class_Type
    private
    integer :: What        ! 1 = IEEE_Quiet_Nan
                           ! 2 = IEEE_Signaling_Nan
  end type IEEE_Class_Type

  type(ieee_class_type), parameter :: ieee_quiet_nan = ieee_class_type(1)
  type(ieee_class_type), parameter :: ieee_signaling_nan = ieee_class_type(2)

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
  ! Private
    
    IEEE_IS_FINITE_S = .FALSE.
    if( isnan(arg) ) RETURN
    if( isinf(arg) ) RETURN
    IEEE_IS_FINITE_S = .TRUE.
  END FUNCTION IEEE_IS_FINITE_S
  
  LOGICAL FUNCTION IEEE_IS_FINITE_D( ARG )
  ! Formal args
    double precision, intent(in) ::          arg
  ! Private
    
    IEEE_IS_FINITE_D = .FALSE.
    if( isnan(arg) ) RETURN
    if( isinf(arg) ) RETURN
    IEEE_IS_FINITE_D = .TRUE.
  END FUNCTION IEEE_IS_FINITE_D
  
  elemental logical function IEEE_Is_NaN_D ( X )
    double precision, intent(in) :: X
    IEEE_Is_NaN_D = isnan(x)
  end function IEEE_Is_NaN_D

  elemental logical function IEEE_Is_NaN_S ( X )
    real, intent(in) :: X
    IEEE_Is_NaN_S = isnan(x)
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

    logical, save :: HaveNaN = .false.
    real, save :: NaN
    character(len=3) :: NaNString = 'NaN'

  ! Private
    if ( .not. haveNaN ) then
      haveNaN = .true.
      read ( NaNString, * ) NaN
    end if

    select case ( class%what )
    case ( ieee_quiet_nan%what )     ! IEEE_Quiet_NaN
      the_value = NaN
    case ( ieee_signaling_nan%what ) ! IEEE_Signaling_NaN
      ! The following is necessary only because I don't know how to generate a
      ! signaling NaN--this obviously fails if any number ever matches this
      the_value = 9.e19
    end select
  END FUNCTION IEEE_VALUE_D

  FUNCTION IEEE_VALUE_S ( X, CLASS ) result ( the_value )
  ! Formal args
    real, intent(in) ::                  X
    type(ieee_class_type), intent(in) :: CLASS
    real ::                              the_value

    logical, save :: HaveNaN = .false.
    real, save :: NaN
    character(len=3) :: NaNString = 'NaN'

  ! Private
    if ( .not. haveNaN ) then
      haveNaN = .true.
      read ( NaNString, * ) NaN
    end if

    select case ( class%what )
    case ( ieee_quiet_nan%what )     ! IEEE_Quiet_NaN
      the_value = NaN
    case ( ieee_signaling_nan%what ) ! IEEE_Signaling_NaN
      ! The following is necessary only because I don't know how to generate a
      ! signaling NaN--this obviously fails if any number ever matches this
      the_value = 9.e19
    end select
  END FUNCTION IEEE_VALUE_S

END MODULE IEEE_ARITHMETIC

!
! $Log$
! Revision 1.4  2004/03/26 21:15:00  pwagner
! FIrst commit
!
! Revision 1.3  2003/07/03 19:26:52  vsnyder
! Remove some comments that are no longer true
!
! Revision 1.2  2003/07/03 19:11:29  vsnyder
! Add IEEE_Support_NaN, IEEE_Is_NaN, IEEE_Class_Type.
! Make IEEE_Quiet_NaN to be of IEEE_Class_Type.
! Make IEEE_Value more like the IEEE TR and Fortran 2003 standard.
!
! Revision 1.1  2001/10/22 18:29:47  pwagner
! First commit
!
