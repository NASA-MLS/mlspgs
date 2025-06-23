! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Parser_Table_m

  ! Type definition for parser tables

  implicit NONE
  private

  public :: Destroy_Parser_Table, Parser_Table_t

  type :: Parser_Table_t
    integer :: IFINAL
    integer :: LSETS
    integer :: NCSETS
    integer :: NLOOKS
    integer :: NSTATE
    integer :: NTERMS
    integer :: NTRANS
    integer :: NUMPRD
    integer :: NVOC
    integer :: TOTAL
    integer :: NBASPROD
    integer, allocatable :: ENT(:)
    integer, allocatable :: FRED(:)
    integer, allocatable :: FTRN(:)
    integer, allocatable :: TRAN(:)
    integer, allocatable :: NSET(:)
    integer, allocatable :: PROD(:)
    integer, allocatable :: LSET(:)
    integer, allocatable :: LS(:)
    integer, allocatable :: LENS(:)
    integer, allocatable :: LHS(:)
    integer, allocatable :: ACT(:)
    integer, allocatable :: VAL(:)
    integer, allocatable :: TEXT(:)
    integer, allocatable :: PROD_IND(:)
    integer, allocatable :: PRODUCTIONS(:)
    integer, allocatable :: INDBAS(:)
    integer, allocatable :: BASPRODS(:)
    integer, allocatable :: DOTS(:)
  end type Parser_Table_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Destroy_Parser_Table ( Intentionally_Not_Used )
    type(parser_table_t), intent(out) :: Intentionally_Not_Used
    ! Intent(OUT) deallocates any allocated allocatable components
  end subroutine Destroy_Parser_Table

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Parser_Table_m

! $Log$
! Revision 2.1  2014/05/21 00:05:49  vsnyder
! Initial commit
!
