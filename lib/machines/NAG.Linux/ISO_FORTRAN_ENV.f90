! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ISO_FORTRAN_ENV

  ! Nonintrinsic version for NAG Fortran v5 for Linux.
  ! It may well work for other versions.

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ISO_FORTRAN_ENV

! $Log$
