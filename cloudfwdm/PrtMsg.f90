
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module PrtMsg

      IMPLICIT NONE
      private
      public :: HEADER

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PrtMsg

! $Log$
! Revision 1.2  2001/09/21 15:51:37  jonathan
! modified F95 version
!
