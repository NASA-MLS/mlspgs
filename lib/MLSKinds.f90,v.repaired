! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSKinds                 ! Kind type parameters for the MLS software
!=============================================================================

  implicit none
  public

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! i1, i2, i4    integer types
! r4, r8        floating point types
! ip, rp        integer, floating point types used in forward model
! rv            floating point type used in vector quantity values
! rm            floating point type used in matrix values

! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Firstly, these are standard numerical kinds, copied from HCP
  ! (again with my change in case, sorry Hugh!)

  integer, public, parameter:: i1=selected_int_kind(2)
  integer, public, parameter:: i2=selected_int_kind(4)
  integer, public, parameter:: i4=selected_int_kind(7)
  integer, public, parameter:: r4=selected_real_kind(5)
  integer, public, parameter:: r8=selected_real_kind(13)

  ! Now choose the precision we want by preference (may automate this through
  ! make later on, with perl or m4 or something).
  ! These are used according to the final letter in the two-letter name:
  ! Final Letter          Context           Suggested Value
  !    ----               -------           ---------------
  !     m                 Matrix            r8  (r4 to save memory)
  !     p                 Forward Model     r8
  !     v                 Vector            r8
  !     t                 Toposet           r4

  integer, public, parameter:: rm=r4
  integer, public, parameter:: rp=r8
  integer, public, parameter:: ip=i4
  integer, public, parameter:: rv=r8
  integer, public, parameter:: rt=r8

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

end module MLSKinds
!=============================================================================

!
! $Log$
! Revision 2.3  2011/08/26 00:21:03  pwagner
! Added numeric type for Interval, Toposet types
!
! Revision 2.2  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2005/10/19 22:52:47  vsnyder
! Initial commit -- move kinds here from MLSCommon
!
