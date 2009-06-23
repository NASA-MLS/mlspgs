
! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
 
module PrtMsg

      IMPLICIT NONE
      private
      public :: HEADER

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains 

      SUBROUTINE HEADER (H)

!=====================================================================
!     J.JIANG, MARCH 8, 2000
!=====================================================================

      INTEGER :: H
!---------------------------------------------------------------------

      IF(H .EQ. 1) THEN
         PRINT*,' '
         PRINT*,'============================================'
         PRINT*,'  >>>>>> START CLOUD FORWARD MODEL <<<<<<<  '
         PRINT*,'============================================'
         PRINT*,' '
         PRINT*,'INPUT MODEL ATMOSPHERE'
      ELSE IF (H .EQ. 2) THEN
         PRINT*,' '
         PRINT*,'COMPUTE CLEAR-SKY EMISSIONS'
      ELSE IF (H .EQ. 3) THEN
         PRINT*,' '
         PRINT*,'COMPUTE CLOUD SCATTERING'
      ELSE IF (H .EQ. 4) THEN
         PRINT*,' '
         PRINT*,'START RADIATIVE TRANSFER CALCULATION...'
      ELSE
         PRINT*,' '
         PRINT*,'COMPLETE !'
      ENDIF

      END SUBROUTINE HEADER

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PrtMsg

! $Log$
! Revision 1.4  2005/06/22 18:27:38  pwagner
! Cant have access declared outside module scope
!
! Revision 1.3  2002/10/08 17:08:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
