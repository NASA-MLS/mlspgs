      SUBROUTINE CLD_TYPE(ITYPE,WC,CHK_CLD,L,LWC1,IWC1,LWC2,IWC2,
     >                                       LWC3,IWC3,LWC4,IWC4)

      INTEGER ITYPE,L
      REAL IWC1(L+1),IWC2(L+1),IWC3(L+1),IWC4(L+1),
     >     LWC1(L+1),LWC2(L+1),LWC3(L+1),LWC4(L+1),
     >     CHK_CLD(L+1)

      REAL WC(2,L+1)     ! CLOUD WATER CONTENT (g/m3)
C---------------------------------------------------------------------

      IF(ITYPE .EQ. 1) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: DEEP-CONVECTIVE SYSTEM '
      ELSE IF (ITYPE .EQ. 2) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: FRONTAL SYSTEM '
      ELSE IF (ITYPE .EQ. 3) THEN
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: ANVILS'
      ELSE 
         PRINT*,' '
         PRINT*,' -> CLOUD-TYPE: THIN-LAYER CIRRUS'
      ENDIF

      DO I=1,L+1

         IF (ITYPE.EQ.1) THEN          ! DEEP-CONVECTIVE  (THICK-CIRRUS) 
            WC(2,I)=0.*LWC1(I)
            WC(1,I)=IWC1(I)
         ELSE IF (ITYPE.EQ.2) THEN     ! FRONTAL SYSTEM   (MID-LAT CIRRUS)
            WC(2,I)=0.*LWC2(I)
            WC(1,I)=IWC2(I)
         ELSE IF (ITYPE.EQ.3) THEN     ! ANVILS           (THIN-CIRRUS1)
            WC(2,I)=0.*LWC3(I)
            WC(1,I)=IWC3(I)
         ELSE IF (ITYPE.EQ.4) THEN     ! THIN-LAYER CIRRUS(NO WATER CLOUDS) 
            WC(2,I)=0.*LWC4(I)
            WC(1,I)=IWC4(I)
         ENDIF
         
         CHK_CLD(I)=WC(1,I)+WC(2,I)
      ENDDO

      RETURN
      END

