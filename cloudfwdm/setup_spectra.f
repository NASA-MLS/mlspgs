      SUBROUTINE SETUP_SPECTRA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,
     >                         GAMMA,N2,MOL,NMOL,NCNT)

      IMPLICIT NONE
      INCLUDE 'spectra.h'
      CHARACTER*80 TITLE
      REAL QQ(3)
      INTEGER NV                          !NO OF SPECTRAL LINES
      INTEGER I, ICNT

      OPEN(31,FILE='spectra.dat', FORM='FORMATTED',STATUS='OLD')
      NMOL=0
      ICNT=0

c      DO WHILE(1)

      DO WHILE(.TRUE.)

         READ(31,4,ERR=100)TITLE
         NMOL=NMOL+1
         NCNT(NMOL)=0
         READ(31,4)TITLE
         READ(31,2)MOL,(QQ(I),I=1,3)
         READ(31,3)NV

         NCNT(NMOL)=NV
         DO I=1,NV
            QLG(1,ICNT+I)=QQ(1)
            QLG(2,ICNT+I)=QQ(2)
            QLG(3,ICNT+I)=QQ(3)
         ENDDO

         READ(31,4)TITLE
         DO I=1,NV
            ICNT=ICNT+1
            READ(31,1)V0(ICNT),GSE(ICNT),IST(ICNT),WTH(ICNT),NTH(ICNT),
     +                DELTA(ICNT),N1(ICNT),GAMMA(ICNT),N2(ICNT)
         ENDDO

      ENDDO

 100  CONTINUE
      CLOSE(31)

 1    FORMAT(1X,F12.4,1X,F10.4,1X,F8.4,1X,F6.3,1X,F5.2,1X,E10.3,1X,
     +       F4.2,1X,E10.3,1X,F5.2)
 2    FORMAT(18X,I4,3X,3(F10.4))
 3    FORMAT(1X,I5)
 4    FORMAT(1X,A80)

      RETURN
      END

