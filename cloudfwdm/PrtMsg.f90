
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

end module PrtMsg

! $Log: PrtMsg.f,v      
