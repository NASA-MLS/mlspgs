
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module CloudProfile

! -------------------------------------------------------------------------  
! MODULE TO GENERATE ICE WATER CONTENT PROFILE IF NEEDED
! -------------------------------------------------------------------------

  use MLSCommon,only: NameLen,    FileNameLen, r8

      IMPLICIT NONE
      private
      public :: CLOUD_MODEL

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------

contains 

      SUBROUTINE CLOUD_MODEL(ITYPE,CHT,YZ,NH,WC)

!=========================================================================C
!  DEFINE VERTICAL PROFILES OF CLOUD ICE-WATER-CONTENT                    C
!  J.JIANG -05/18/2001                                                    C
!          -10/05/2001, MODIFIED TO FIT CLOUD RETREVIAL REQUIREMENTS      C
!=========================================================================C

      CHARACTER :: ITYPE

      INTEGER :: NH,I

      REAL(r8) :: YZ(NH)
      REAL(r8) :: WC(2,NH)

      REAL(r8) :: CLD_TOP
      REAL(r8) :: CLD_BASE
      REAL(r8) :: UPPER_LAG
      REAL(r8) :: LOWER_LAG
      REAL(r8) :: HT, CHT, WCscale

!--------------------------------------------------------------------

      IF (CHT .EQ. 0._r8) THEN
         HT = 16000. + CHT*1000.
      ELSE
         HT=CHT
      ENDIF

      WCscale=0.1

      IF(ITYPE .EQ. 'C') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: DEEP-CONVECTIVE SYSTEM '

         CLD_TOP   = HT + 1000.         ! CONVECTIVE CLOUD-TOP
         CLD_BASE  = HT - 11600.        ! CONVECTIVE CLOUD-BASE
         UPPER_LAG = HT - 8000.
         LOWER_LAG = HT - 10000.

         DO I=1,NH

            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/5000.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF

            ! LIQUID WATER CLOUD
            
            IF (YZ(I).GT.5000. .AND. YZ(I).LT. 10000.) THEN
               WC(2,I) = 0.1*0.
            ENDIF
         ENDDO


      ELSE IF (ITYPE .EQ. 'F') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: FRONTAL SYSTEM '

         CLD_TOP   = HT - 5000.         ! FRONTAL CLOUD-TOP
         CLD_BASE  = HT - 11000.        ! FRONTAL CLOUD-BASE
         UPPER_LAG = HT - 7000.
         LOWER_LAG = HT - 9000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 'A') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: ANVILS'

         CLD_TOP   = HT - 4000.         !ANVIL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !ANVIL CLOUD-BASE
         UPPER_LAG = HT - 5000.
         LOWER_LAG = HT - 6000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/600.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 'T') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: THIN-LAYER CIRRUS'

         CLD_TOP   = HT + 500.         !1KM THICK CLOUD LAYER
         CLD_BASE  = HT - 500.
         UPPER_LAG = HT + 200.
         LOWER_LAG = HT - 200.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-HT)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(HT-YZ(I))/500.)
            ENDIF
            WC(2,I)=0.0
         ENDDO
 
      ELSE

        PRINT*, 'NO CLOUDS'
       
        DO I=1,NH
            WC(1,I)=0.0
            WC(2,I)=0.0
        ENDDO

      ENDIF


        DO I=1,NH
            WC(1,I)=WCscale*WC(1,I)
            WC(2,I)=WCscale*WC(2,I)
        ENDDO


!------------------------------------------------------

      END SUBROUTINE CLOUD_MODEL

end module CloudProfile

! $Log$     



