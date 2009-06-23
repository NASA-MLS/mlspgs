! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ISO_FORTRAN_ENV

  ! Nonintrinsic version for Lahey/Fujitsu Fortran v6.10c for Linux.
  ! It may well work for other versions.

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, parameter :: Character_Storage_Size = 8
  integer, parameter :: Error_Unit = 0
  integer, parameter :: File_Storage_Size = 8
  integer, parameter :: Input_Unit = 5
  integer, parameter :: IOSTAT_END = -1
  integer, parameter :: IOSTAT_EOR = -2
  integer, parameter :: Numeric_Storage_Size = 32
  integer, parameter :: Output_Unit = 6

contains

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ISO_FORTRAN_ENV

! $Log$
! Revision 1.2  2005/06/22 20:26:22  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.1  2004/10/06 00:20:58  vsnyder
! Initial commit
!
