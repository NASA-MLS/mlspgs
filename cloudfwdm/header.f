      SUBROUTINE HEADER (H)

      INTEGER H
C---------------------------------------------------------------------

      IF(H .EQ. 1) THEN
         PRINT*,' '
         PRINT*,'============================================'
         PRINT*,'EOS MLS MICROWAVE RADIATIVE TRANSFER PROGRAM'
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

      RETURN
      END

