! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_NZ_m

  ! Dump the NZ and NNZ arrays from Get_Eta_Matrix.
  ! This is not in Get_Eta_Matrix so as to avoid creating a dependence
  ! therein upon the Dump_0 and Output_m modules.

  implicit NONE

  private

  public :: Dump_NZ

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Dump_NZ ( NZ, NNZ, What )

    use Dump_0, only: Dump
    use Output_m, only: Output, NewLine

    integer, intent(in) :: NZ(:,:) ! Nonzero locations in Eta
    integer, intent(in) :: NNZ(:)  ! Number of nonzeroes in each column
    character(len=*), intent(in), optional :: What

    integer :: I

    call output ( 'Dump of nonzero structure' )
    if ( present(what) ) call output ( trim(what) )
    call newLine

    do i = 1, size(nnz)
      call output ( i, before='Nonzeros in column ', advance='yes' )
      call dump ( nz(:nnz(i),i) )
    end do

  end subroutine Dump_NZ

!=========================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_NZ_m

! $Log$
! Revision 2.2  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2007/09/07 01:34:00  vsnyder
! Initial commit
!
