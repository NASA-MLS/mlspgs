! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject bto U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module To_Log_Basis_m
!=============================================================================

! Copy an array from linear to log basis.

  implicit NONE
  private

  ! Public procedures:

  public :: To_Log_Basis

  interface To_Log_Basis
    module procedure To_Log_Basis_0,  To_Log_Basis_1, To_Log_Basis_2, &
      & To_Log_Basis_3
  end interface

  !------------- RCS Ident Info (more below in not_used_here) ----------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  elemental subroutine To_Log_Basis_0 ( Value, Min_Val )
    use VectorsModule, only: RV
    real(rv), intent(inout) :: Value ! Converted to log(linear) on output
    real(rv), intent(in) :: Min_Val
    if ( value < min_val ) then
      value = log(min_val)
    else
      value = log(value)
    end if
  end subroutine To_Log_Basis_0

  pure &
  subroutine To_Log_Basis_1 ( Value, Min_Val )
    use VectorsModule, only: RV
    real(rv), intent(inout) :: Value(:)
    real(rv), intent(in) :: Min_Val
    where ( value < min_val )
      value = log(min_val)
    else where
      value = log(value)
    end where
  end subroutine To_Log_Basis_1

  pure &
  subroutine To_Log_Basis_2 ( Value, Min_Val )
    use VectorsModule, only: RV
    real(rv), intent(inout) :: Value(:,:)
    real(rv), intent(in) :: Min_Val
    where ( value < min_val )
      value = log(min_val)
    else where
      value = log(value)
    end where
  end subroutine To_Log_Basis_2

  pure &
  subroutine To_Log_Basis_3 ( Value, Min_Val )
    use VectorsModule, only: RV
    real(rv), intent(inout) :: Value(:,:,:)
    real(rv), intent(in) :: Min_Val
    where ( value < min_val )
      value = log(min_val)
    else where
      value = log(value)
    end where
  end subroutine To_Log_Basis_3

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module To_Log_Basis_m

! $Log$
! Revision 2.1  2017/01/26 00:51:26  vsnyder
! Initial commit
!
