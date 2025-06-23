! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module SDPToolkit     ! Substitute for the essential toolkit routines
!=============================================================================

  implicit none

  public::PGS_SMF_GenerateStatusReport
  !---------------------------- RCS Ident Info -------------------------------
  character (len = 256),parameter,private :: Id = &
       "$Id$"
  character ( len = *), parameter,private :: &
       ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------
  ! This module gives substitutes for the essential toolkit routines used by
  ! library code used in both the toolkit and non toolkit environment.
  ! This is HCPs personal version, which has the same module name and  
  ! filename as the true SDP Toolkit module. I think this deserves a 
  ! separate directory from mlspgs/lib, so one can tell the compiler to look
  ! there first for source files

contains

  function PGS_SMF_GenerateStatusReport(message) result (GenerateStatusReport)
    character (len=*), intent(in) :: message
    integer :: GenerateStatusReport

    print*,message
    GenerateStatusReport=0
  end function PGS_SMF_GenerateStatusReport

!=============================================================================
end module SDPToolkit
!=============================================================================

