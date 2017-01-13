! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_Row_Sparse_m

  implicit NONE
  private

  ! Dump a two-dimensional array that is row sparse, printing only the nonzeros.
  ! If it actually has more dimensions packed into the second dimension, the
  ! Bounds2 argument gives their upper bounds (lower bounds assumed to be 1),
  ! which are then used to calculate those subscripts instead of simply showing
  ! a column number.

  public :: Dump_Row_Sparse

  interface Dump_Row_Sparse
    module procedure Dump_Row_Sparse_D, Dump_Row_Sparse_S
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Dump_Row_Sparse_d ( Array, Name, Width, Format, Bounds2 )
    use Array_Stuff, only: Subscripts
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: Blanks, NewLine, Output
    double precision, intent(in) :: Array(:,:)
    character(*), intent(in), optional :: Name
    integer, intent(in), optional :: Width     ! Array elements per line
    character(*), intent(in), optional :: Format
    integer, intent(in), optional :: Bounds2(:) ! Bounds of 2nd dim array
    include "Dump_Row_Sparse.f9h"
  end subroutine Dump_Row_Sparse_d

  subroutine Dump_Row_Sparse_s ( Array, Name, Width, Format, Bounds2 )
    use Array_Stuff, only: Subscripts
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: Blanks, NewLine, Output
    real, intent(in) :: Array(:,:)
    character(*), intent(in), optional :: Name
    integer, intent(in), optional :: Width     ! Array elements per line
    character(*), intent(in), optional :: Format
    integer, intent(in), optional :: Bounds2(:) ! Bounds of 2nd dim array
    include "Dump_Row_Sparse.f9h"
  end subroutine Dump_Row_Sparse_s

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_Row_Sparse_m

! $Log$
! Revision 2.1  2017/01/13 21:01:58  vsnyder
! Initial commit
!
