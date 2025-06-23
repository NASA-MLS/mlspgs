
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module HT_M ! From Lawson and Hanson, converted to Fortran 90
!=============================================================================

  implicit NONE
  private

  public :: HTGEN, SHTGENM, SHTVV, DHTGENM, DHTVV

  interface HTGEN
    module procedure SHTGENM, SHTVV, DHTGENM, DHTVV
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine SHTGENM ( MODE, LPIVOT, L1, M, U, UPARAM, C, COLC )

    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in)     :: MODE, LPIVOT, L1, M
    real(rk), intent(inout) :: U(:), UPARAM, C(:,:)
    logical, intent(in)     :: COLC
    include 'htgenm.f9h'

  end subroutine SHTGENM

  subroutine SHTVV ( MODE, LPIVOT, L1, M, U, UPARAM, C )

    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0e0)
    integer, intent(in)     :: MODE, LPIVOT, L1, M
    real(rk), intent(inout) :: U(:), UPARAM, C(:)
    include 'htvv.f9h'

  end subroutine SHTVV

  subroutine DHTGENM ( MODE, LPIVOT, L1, M, U, UPARAM, C, COLC )

    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in)     :: MODE, LPIVOT, L1, M
    real(rk), intent(inout) :: U(:), UPARAM, C(:,:)
    logical, intent(in)     :: COLC
    include 'htgenm.f9h'

  end subroutine DHTGENM

  subroutine DHTVV ( MODE, LPIVOT, L1, M, U, UPARAM, C )

    use NRM2_M, only: NRM2

    integer, parameter :: RK = kind(0.0d0)
    integer, intent(in)     :: MODE, LPIVOT, L1, M
    real(rk), intent(inout) :: U(:), UPARAM, C(:)
    include 'htvv.f9h'

  end subroutine DHTVV

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HT_M

! $Log$
! Revision 2.2  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2006/03/22 02:07:04  vsnyder
! Initial commit
!
