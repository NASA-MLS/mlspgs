! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_IEEE_NaN_m

  use, intrinsic :: IEEE_Arithmetic, only: IEEE_Signaling_NaN, IEEE_Value, &
    IEEE_Support_DataType
  implicit NONE
  private

  public :: Fill_IEEE_NaN, Get_IEEE_NaN

  ! Fill up to 26 arguments, all the same rank, with NaN
  ! It might seem to be nice to use an elemental subroutine, but
  ! that would require all the arguments to be conforming, i.e., all
  ! either scalars, or arrays of the same shape.
  interface Fill_IEEE_NaN
    module procedure Fill_IEEE_NaN_D_0D, Fill_IEEE_NaN_S_0D
    ! The next eight are needed because Intel ifort 14.0.0.080
    ! refused to compile the elemental ones.
    module procedure Fill_IEEE_NaN_D_1D, Fill_IEEE_NaN_S_1D
    module procedure Fill_IEEE_NaN_D_2D, Fill_IEEE_NaN_S_2D
    module procedure Fill_IEEE_NaN_D_3D, Fill_IEEE_NaN_S_3D
    module procedure Fill_IEEE_NaN_D_4D, Fill_IEEE_NaN_S_4D
  end interface

  interface Get_IEEE_NaN
    module procedure Get_IEEE_NaN_D, Get_IEEE_NaN_S
  end interface

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains

! Intel ifort 14.0.0.080 refused to compile these, claiming that
! IEEE_Support_DataType is not pure.
!   elemental subroutine Fill_IEEE_NaN_D ( X )
!     double precision, intent(out) :: X
!     X = merge ( IEEE_Value(x,IEEE_Signaling_NaN), 0.0d0, &
!               & IEEE_Support_DataType(x) )
!   end subroutine Fill_IEEE_NaN_D
! 
!   elemental subroutine Fill_IEEE_NaN_S ( X )
!     real, intent(out) :: X
!     X = merge ( IEEE_Value(x,IEEE_Signaling_NaN), 0.0e0, &
!               & IEEE_Support_DataType(x) )
!   end subroutine Fill_IEEE_NaN_S

  subroutine Fill_IEEE_NaN_D_0D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_D_0D

  subroutine Fill_IEEE_NaN_D_1D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), dimension(:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_D_1D

  subroutine Fill_IEEE_NaN_D_2D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), dimension(:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_D_2D

  subroutine Fill_IEEE_NaN_D_3D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), dimension(:,:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_D_3D

  subroutine Fill_IEEE_NaN_D_4D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), dimension(:,:,:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_D_4D

  subroutine Fill_IEEE_NaN_S_0D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_S_0D

  subroutine Fill_IEEE_NaN_S_1D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), dimension(:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_S_1D

  subroutine Fill_IEEE_NaN_S_2D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), dimension(:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_S_2D

  subroutine Fill_IEEE_NaN_S_3D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), dimension(:,:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_S_3D

  subroutine Fill_IEEE_NaN_S_4D ( A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), dimension(:,:,:,:), intent(out) :: &
      & A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    optional :: &
      & B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z
    include 'Fill_NaN.f9h'
  end subroutine Fill_IEEE_NaN_S_4D

  double precision function Get_IEEE_NaN_D ( X )
    double precision, intent(in) :: X
    Get_IEEE_NaN_D = merge ( IEEE_Value(x,IEEE_Signaling_NaN), 0.0d0, &
                           & IEEE_Support_DataType(x) )
  end function Get_IEEE_NaN_D

  real function Get_IEEE_NaN_S ( X )
    real, intent(in) :: X
    Get_IEEE_NaN_S = merge ( IEEE_Value(x,IEEE_Signaling_NaN), 0.0e0, &
                           & IEEE_Support_DataType(x) )
  end function Get_IEEE_NaN_S

  ! ------------------------------------------------  Not_Used_Here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_IEEE_NaN_m

! $Log$
! Revision 2.1  2014/07/04 00:10:08  vsnyder
! Initial commit
!
