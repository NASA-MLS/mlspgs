! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE IEEE_ARITHMETIC              ! Common utilities for the MLSL1 program
!=============================================================================

!  use SUNSOWNIEEE, only: ir_isnan, ir_isinf, r_quiet_nan

  implicit NONE
  private

  ! These are the ieee functions and constants needed in MLS
  ! that Sun fails to supply.

  public :: IEEE_Class_Type
  public :: IEEE_Is_finite, IEEE_Is_NaN, IEEE_Quiet_NaN, IEEE_Support_NaN
  public :: IEEE_Value

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains glue routines for Sun's own f95 compiler
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

  elemental logical function IEEE_IS_FINITE_S( ARG )
  ! Formal args
    real, intent(in) ::          arg
  ! Private
    
    IEEE_IS_FINITE_S = ( abs(arg) <= Huge(arg) )
  END FUNCTION IEEE_IS_FINITE_S
  
  elemental logical function IEEE_IS_FINITE_D( ARG )
  ! Formal args
    double precision, intent(in) ::          arg
  ! Private
    
    IEEE_IS_FINITE_D = ( abs(arg) <= Huge(arg) )
  END FUNCTION IEEE_IS_FINITE_D
  
  elemental logical function IEEE_Is_Inf_io_D ( X ) result(res)
    double precision, intent(in) :: X
    character(len=80) :: reschar
    write(reschar, *) x
    res = ( index(reschar, 'Inf') > 0 .or. &
      & index(reschar, 'inf') > 0 .or. index(reschar, 'INF') > 0 )
  end function IEEE_Is_Inf_io_D

  elemental logical function IEEE_Is_Inf_io_S ( X ) result(res)
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

  ! The following is made necessary only because Sun's own f95
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

  ! The following is made necessary only because Sun's own f95
  ! fails to comply--this obviously fails if any number ever matches this
    select case ( class%what )
    case ( ieee_quiet_nan%what )     ! IEEE_Quiet_NaN
      the_value = -9.e19
    case ( ieee_signaling_nan%what ) ! IEEE_Signaling_NaN
      the_value = 9.e19
    end select
  END FUNCTION IEEE_VALUE_S

! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE IEEE_ARITHMETIC

!
! $Log$
! Revision 1.9  2009/06/23 19:58:53  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 1.8  2005/06/22 20:27:58  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.7  2005/05/12 20:39:11  pwagner
! Made ieee_is_finite elemental
!
! Revision 1.6  2004/05/24 18:21:55  pwagner
! IEEE_Is_NaN can handle both single, double precision
!
! Revision 1.5  2003/07/03 19:25:24  vsnyder
! Remove some comments that are no longer true
!
! Revision 1.4  2003/07/03 19:11:57  vsnyder
! Add IEEE_Support_NaN, IEEE_Is_NaN, IEEE_Class_Type.
! Make IEEE_Quiet_NaN to be of IEEE_Class_Type.
! Make IEEE_Value more like the IEEE TR and Fortran 2003 standard.
!
! Revision 1.3  2001/09/18 17:18:12  pwagner
! Strictly so the darn thing compiles--let others (or Sun) make it work
!
! Revision 1.2  2001/09/14 17:35:33  pwagner
! Fixed an obvious blunder; expanded some comments
!
! Revision 1.1  2001/09/13 23:18:37  pwagner
! First commit
!
