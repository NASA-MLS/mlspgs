! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Functions

! Provide parameters to identify functions.
! Put names of functions into the declaration table, with their
! corresponding parameter in the "units" field.

! We may eventually want to put a tree for each function into init_tables
! to specify the allowed argument types and the result type, and put some
! checking in tree_checker.

  implicit NONE
  public

  integer, parameter :: F_Cholesky = 1
  integer, parameter :: F_Transpose = F_Cholesky + 1

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Init_Functions

    call declare_func ( 'cholesky', f_cholesky )
    call declare_func ( 'transpose', f_transpose )

  contains
    subroutine Declare_Func ( String, Index )
    
      use DECLARATION_TABLE, only: DECLARE, FUNCTION
      use INTRINSIC, only: ADD_IDENT
      use TREE, only: NULL_TREE

      character(len=*), intent(in) :: String ! The text of the function
      integer, intent(in) :: Index           ! The index number for the function
      call declare ( add_ident(string), 0.0d0, function, index, null_tree )
    end subroutine Declare_Func
  end subroutine Init_Functions

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Functions

! $Log$
! Revision 2.1  2004/01/17 03:05:06  vsnyder
! Initial commit
!
