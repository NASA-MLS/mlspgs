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
    module procedure Remap_3d_Char, Remap_3d_Double, Remap_3d_Real
    module procedure Remap_4d_Char, Remap_4d_Double, Remap_4d_Real
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

! include lines aren't used here because Intel ifort 12.0-1230 didn't want
! to process the #ifdef directives that would be in them.

  subroutine Remap_2d_Char ( A, B, Shape, Lbounds )
    character, pointer :: A(:), B(:,:)
    integer, intent(in) :: Shape(2)
    integer, intent(in), optional :: Lbounds(2)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2)) => a
    else
      b(1:shape(1),1:shape(2)) => a
    end if
#endif

  end subroutine Remap_2d_Char

  subroutine Remap_2d_Double ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:)
    integer, intent(in) :: Shape(2)
    integer, intent(in), optional :: Lbounds(2)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2)) => a
    else
      b(1:shape(1),1:shape(2)) => a
    end if
#endif

  end subroutine Remap_2d_Double

  subroutine Remap_2d_Real ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:)
    integer, intent(in) :: Shape(2)
    integer, intent(in), optional :: Lbounds(2)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2)) => a
    else
      b(1:shape(1),1:shape(2)) => a
    end if
#endif

  end subroutine Remap_2d_Real

  subroutine Remap_3d_Char ( A, B, Shape, Lbounds )
    character, pointer :: A(:), B(:,:,:)
    integer, intent(in) :: Shape(3)
    integer, intent(in), optional :: Lbounds(3)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3)) => a
    end if
#endif

  end subroutine Remap_3d_Char

  subroutine Remap_3d_Double ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:,:)
    integer, intent(in) :: Shape(3)
    integer, intent(in), optional :: Lbounds(3)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3)) => a
    end if
#endif

  end subroutine Remap_3d_Double

  subroutine Remap_3d_Real ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:,:)
    integer, intent(in) :: Shape(3)
    integer, intent(in), optional :: Lbounds(3)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3)) => a
    end if
#endif

  end subroutine Remap_3d_Real

  subroutine Remap_4d_Char ( A, B, Shape, Lbounds )
    character, pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: Shape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3),lbounds(4):shape(4)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3),1:shape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Char

  subroutine Remap_4d_Double ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: Shape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3),lbounds(4):shape(4)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3),1:shape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Double

  subroutine Remap_4d_Real ( A, B, Shape, Lbounds )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), pointer :: A(:), B(:,:,:,:)
    integer, intent(in) :: Shape(4)
    integer, intent(in), optional :: Lbounds(4)
#ifdef CMAP
    type(c_ptr) :: C
    c = c_loc(a(1))
    call c_f_pointer ( c, b, shape )
#ifdef LBOUND
    if ( present(lbounds) ) b(lbounds(1):,lbounds(2):,lbounds(3):,lbounds(4):) => b
#else
    if ( present(lbounds) ) call MLSMessage ( MLSMSG_Crash, &
      & "Pointer low bounds setting not supported by compiler version", &
      & moduleName )
#endif
#else
    if ( present(lbounds) ) then
      b(lbounds(1):shape(1),lbounds(2):shape(2),lbounds(3):shape(3),lbounds(4):shape(4)) => a
    else
      b(1:shape(1),1:shape(2),1:shape(3),1:shape(4)) => a
    end if
#endif

  end subroutine Remap_4d_Real

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
! Revision 2.2  2012/07/10 00:05:21  vsnyder
! More #ifdef stuff to handle rank remapping using C instead of Fortran
! 2003.  Add LBOUND to turn on setting lower bounds in the C case.  Assume
! setting lower bounds works if Fortran remapping works.
!
! Revision 2.1  2012/07/07 01:59:41  vsnyder
! Initial commit
!
