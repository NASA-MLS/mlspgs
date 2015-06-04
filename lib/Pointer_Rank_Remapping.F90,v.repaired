! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Pointer_Rank_Remapping

! This module contains subroutines to re-map the rank of a 1-D pointer to
! higher rank.

#ifdef CMAP
  use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
#endif
  use MLSMessageModule, only: MLSMSG_Crash, MLSMessage
  implicit NONE
  private

  public :: Remap

  interface Remap
    module procedure Remap_2d_Char, Remap_2d_Double, Remap_2d_Real
    module procedure Remap_2d_to_1d_Double, Remap_2d_to_1d_Real
    module procedure Remap_2d_to_3d_Double, Remap_2d_to_3d_Real
    module procedure Remap_2d_to_4d_Double, Remap_2d_to_4d_Real
    module procedure Remap_3d_Char, Remap_3d_Double, Remap_3d_Real
    module procedure Remap_3d_to_2d_Double, Remap_3d_to_2d_Real
    module procedure Remap_4d_Char, Remap_4d_Double, Remap_4d_Real
  ! For an allocatable, target argument
    module procedure RemapA_2d_Char, RemapA_2d_Double, RemapA_2d_Real
    module procedure RemapA_3d_Char, RemapA_3d_Double, RemapA_3d_Real
    module procedure RemapA_4d_Char, RemapA_4d_Double, RemapA_4d_Real
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

! include lines aren't used here because Intel ifort 12.0-1230 didn't want
! to process the #ifdef directives that would be in them.

  subroutine Remap_2d_Char ( A, B, TheShape, Lbounds )
    character, pointer :: A(:), B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
integer :: One, Two ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
!      b(1:theShape(1),1:theShape(2)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
b(1:one,1:two) => a
    end if
#endif

  end subroutine Remap_2d_Char

  subroutine Remap_2d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
integer :: One, Two ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
!      b(1:theShape(1),1:theShape(2)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
b(1:one,1:two) => a
    end if
#endif

  end subroutine Remap_2d_Double

  subroutine Remap_2d_to_1d_Double ( A, B, TheShape, Lbound )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:,:), B(:)
    integer, intent(in) :: TheShape
    integer, intent(in), optional :: Lbound
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, [ theShape ] )
  ! Eventually, do this:
  ! b(1:theShape) => a
    if ( present(lbound) ) b(lbound:) => b
  end subroutine Remap_2d_to_1d_Double

  subroutine Remap_2d_to_3d_Double ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:,:), B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
  end subroutine Remap_2d_to_3d_Double

  subroutine Remap_2d_to_4d_Double ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:,:), B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
  end subroutine Remap_2d_to_4d_Double

  subroutine Remap_2d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
integer :: One, Two ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
!      b(1:theShape(1),1:theShape(2)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
b(1:one,1:two) => a
    end if
#endif

  end subroutine Remap_2d_Real

  subroutine Remap_2d_to_1d_Real ( A, B, TheShape, Lbound )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:,:), B(:)
    integer, intent(in) :: TheShape
    integer, intent(in), optional :: Lbound
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, [ theShape ] )
  ! Eventually, do this:
  ! b(1:size(a)) => a
    if ( present(lbound) ) b(lbound:) => b
  end subroutine Remap_2d_to_1d_Real

  subroutine Remap_2d_to_3d_Real ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:,:), B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
  end subroutine Remap_2d_to_3d_Real

  subroutine Remap_2d_to_4d_Real ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:,:), B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
  end subroutine Remap_2d_to_4d_Real

  subroutine Remap_3d_Char ( A, B, TheShape, Lbounds )
    character, pointer :: A(:), B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine Remap_3d_Char

  subroutine Remap_3d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine Remap_3d_Double

  subroutine Remap_3d_to_2d_Double ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:,:,:), B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
  end subroutine Remap_3d_to_2d_Double

  subroutine Remap_3d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine Remap_3d_Real

  subroutine Remap_3d_to_2d_Real ( A, B, TheShape, Lbounds )
    use, intrinsic :: ISO_C_Binding, only: C_F_Pointer, C_LOC, C_PTR
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:,:,:), B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
    type(c_ptr) :: C
    ! For now, use the kludge of C_Loc and C_F_Pointer because some
    ! compilers don't allow targets with rank greater than one
    c = c_loc(a(1,1,1))
    call c_f_pointer ( c, b, theShape )
  ! Eventually, do this:
  ! b(1:theShape(1),1:theShape(2)) => a
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
  end subroutine Remap_3d_to_2d_Real

  subroutine Remap_4d_Char ( A, B, TheShape, Lbounds )
    character, pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Char

  subroutine Remap_4d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Double

  subroutine Remap_4d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Real

  subroutine RemapA_2d_Char ( A, B, TheShape, Lbounds )
    character, allocatable, target :: A(:)
    character, pointer :: B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
integer :: One, Two ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
!      b(1:theShape(1),1:theShape(2)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
b(1:one,1:two) => a
    end if
#endif

  end subroutine RemapA_2d_Char

  subroutine RemapA_2d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
      b(1:theShape(1),1:theShape(2)) => a
    end if

  end subroutine RemapA_2d_Double

  subroutine RemapA_2d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:)
    integer, intent(in) :: TheShape(2)
    integer, intent(in), optional :: Lbounds(2)
integer :: One, Two ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2)) => a
    else
!      b(1:theShape(1),1:theShape(2)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
b(1:one,1:two) => a
    end if
#endif

  end subroutine RemapA_2d_Real

  subroutine RemapA_3d_Char ( A, B, TheShape, Lbounds )
    character, allocatable, target :: A(:)
    character, pointer :: B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine RemapA_3d_Char

  subroutine RemapA_3d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine RemapA_3d_Double

  subroutine RemapA_3d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:,:)
    integer, intent(in) :: TheShape(3)
    integer, intent(in), optional :: Lbounds(3)
integer :: One, Two, Three ! A kludge to get around a NAG runtime-check bug
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3)) => a
    else
!      b(1:theShape(1),1:theShape(2),1:theShape(3)) => a
! This is a replacement to get around a bug in NAG 5.2(723)
one = theShape(1)
two = theShape(2)
three = theShape(3)
b(1:one,1:two,1:three) => a
    end if
#endif

  end subroutine RemapA_3d_Real

  subroutine RemapA_4d_Char ( A, B, TheShape, Lbounds )
    character, allocatable, target :: A(:)
    character, pointer :: B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine RemapA_4d_Char

  subroutine RemapA_4d_Double ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine RemapA_4d_Double

  subroutine RemapA_4d_Real ( A, B, TheShape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), allocatable, target :: A(:)
    real(rk), pointer :: B(:,:,:,:)
    integer, intent(in) :: TheShape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, theShape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):theShape(1),lbounds(2):theShape(2),lbounds(3):theShape(3),lbounds(4):theShape(4)) => a
    else
      b(1:theShape(1),1:theShape(2),1:theShape(3),1:theShape(4)) => a
    end if
#endif

  end subroutine RemapA_4d_Real

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Pointer_Rank_Remapping

! $Log$
! Revision 2.5  2015/06/04 01:55:49  vsnyder
! Add remappings from other than 1D
!
! Revision 2.4  2015/06/02 23:54:31  vsnyder
! Provide for remapping allocatable arrays
!
! Revision 2.3  2012/07/19 03:39:20  vsnyder
! Replace 'shape' by 'theShape', work around a bug in NAG 5.2(723)
!
! Revision 2.2  2012/07/10 00:05:21  vsnyder
! More #ifdef stuff to handle rank remapping using C instead of Fortran
! 2003.  Add LBOUND to turn on setting lower bounds in the C case.  Assume
! setting lower bounds works if Fortran remapping works.
!
! Revision 2.1  2012/07/07 01:59:41  vsnyder
! Initial commit
!
