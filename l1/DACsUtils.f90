! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE DACsUtils
!=============================================================================

  USE MLSL1Utils, ONLY: BigEndianStr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ExtractDACSdata, UncompressDACSdata, ProcessDACSdata, InitDACS_FFT

  INTEGER :: plan129

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE InitDACS_FFT
!=============================================================================

! FFTW Parameters

    INTEGER, PARAMETER :: FFTW_FORWARD = -1
    INTEGER, PARAMETER :: FFTW_BACKWARD = 1
    INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX = -1
    INTEGER, PARAMETER :: FFTW_COMPLEX_TO_REAL = 1
    INTEGER, PARAMETER :: FFTW_ESTIMATE = 0
    INTEGER, PARAMETER :: FFTW_MEASURE = 1
    INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE = 0
    INTEGER, PARAMETER :: FFTW_IN_PLACE = 8
    INTEGER, PARAMETER :: FFTW_USE_WISDOM = 16

! Set up plan for 129 channel DACS:

    CALL rfftw_f77_create_plan (plan129, 256, FFTW_REAL_TO_COMPLEX, &
         FFTW_MEASURE)

  END SUBROUTINE InitDACS_FFT

!=============================================================================
  SUBROUTINE ExtractDACSdata (rawdata, K, D, TP, DIO, LO)
!=============================================================================

    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: rawdata
    INTEGER, DIMENSION(:), INTENT(OUT) :: K
    INTEGER, INTENT(OUT) :: D(4), TP, DIO, LO

    INTEGER :: i, nchans, offset

    K = 0

    IF (SIZE (rawdata) > 400) THEN   ! contains 129 channels of data
       nchans = 129
    ELSE                             ! truncated DACS data
       nchans = 82
    ENDIF

    DO i = 1, nchans*3, 3            ! every 3 bytes (MSB)
       K(i/3+1) = BigEndianStr (rawdata(i)//rawdata(i+1)//rawdata(i+2))
    ENDDO

    offset = nchans * 3 + 1

    DO i = 1, 4
       D(i) = BigEndianStr (rawdata(offset)//rawdata(offset+1)// &
            rawdata(offset+2))
       offset = offset + 3           ! next 3 bytes
    ENDDO

    TP = BigEndianStr (rawdata(offset)//rawdata(offset+1))

    IF (nchans > 82) THEN
       DIO = BigEndianStr (rawdata(offset+2))
       LO = BigEndianStr (rawdata(offset+3))
    ENDIF

  END SUBROUTINE ExtractDACSdata

!=============================================================================
  SUBROUTINE UncompressDACSdata (rawdata, C, D, TP, DIO, LO, Zlag)
!=============================================================================

    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: rawdata
    INTEGER, INTENT(OUT) :: C(*), D(4), TP, DIO, Zlag, LO

    INTEGER :: i, i12, i8, rindx

!! Masks for extending sign bits:

    INTEGER :: bitmask7, signmask7, bitmask11, signmask11, upper4bits
    DATA bitmask7 / z'80' /
    DATA signmask7 / z'FFFFFF00' /
    DATA bitmask11 / z'800' /
    DATA signmask11 / z'FFFFF000' /
    DATA upper4bits / z'FFF' /

! Get the 12 bit quantities

    DO i = 1, 31, 2      ! odd indexes
       rindx = i + i / 2
       i12 = BigEndianStr (rawdata(rindx)//rawdata(rindx+1))
       i12 = ISHFT (i12, -4)           ! shift 4 bits to the right
       IF (IAND (i12, bitmask11) /= 0) THEN
          i12 = IOR (i12, signmask11)
       ENDIF
       C(i) = i12
    ENDDO

    DO i = 2, 32, 2      ! even indexes
       rindx = i + (i - 1) / 2
       i12 = BigEndianStr (rawdata(rindx)//rawdata(rindx+1))
       i12 = IAND (i12, upper4bits)        ! mask off upper 4 bits
       IF (IAND (i12, bitmask11) /= 0) THEN
          i12 = IOR (i12, signmask11)
       ENDIF
       C(i) = i12
    ENDDO

! Get the 8 bit quantities

    DO i = 33, 128
       rindx = i + 16
       i8 = BigEndianStr (rawdata(rindx))
       IF (IAND (i8, bitmask7) /= 0) THEN
          i8 = IOR (i8, signmask7)
       ENDIF
       C(i) = i8
    ENDDO

! Get the D array values

    DO i = 1, 4
       rindx = 145 + (i - 1) * 3
       D(i) = BigEndianStr (rawdata(rindx)//rawdata(rindx+1)//rawdata(rindx+2))
    ENDDO

    rindx = rindx + 3
    TP = BigEndianStr (rawdata(rindx)//rawdata(rindx+1))

    IF (SIZE (rawdata) > 158) THEN
       DIO = BigEndianStr (rawdata(rindx+2))
       Zlag = BigEndianStr (rawdata(rindx+3))
       LO = BigEndianStr (rawdata(rindx+4))
    ELSE
       DIO = -1
       Zlag = -1
       LO = -1
    ENDIF

  END SUBROUTINE UncompressDACSdata

!=============================================================================
  SUBROUTINE ProcessDACSdata
!=============================================================================

    USE MLSL1Common, ONLY: DACSnum, DACSchans, MaxMIFs
    USE L0_sci_tbls, ONLY: DACS_MAF, SciMAF

    INTEGER :: DACSno, MIFno, nchans
    LOGICAL :: GoodDACS
    REAL :: R2a(DACSchans)

    CALL FixLostCarryBits  ! fix state counters for the MAF

    DO MIFno = 0, (MaxMIFs - 1)
       DO DACSno = 1, DACSnum
          GoodDACS = (SUM (DACS_MAF(MIFno)%D(:,DACSno)) > 0)
          IF (GoodDACS) THEN
             nchans = 129
             IF ((.NOT. DACS_MAF(MIFno)%Compressed) .AND. &
                  (MOD (DACSno, 2) == 0)) nchans = 82
             CALL UnpackDACSdata (DACS_MAF(MIFno)%C_K(:,DACSno), &
                  DACS_MAF(MIFno)%D(:,DACSno), DACS_MAF(MIFno)%Compressed, &
                  nchans, R2a)
             CALL ProcessUnpackedDACS (DACS_MAF(MIFno)%D(:,DACSno), R2a, &
                  nchans, DACS_MAF(MIFno)%TP(DACSno), &
                  SciMAF(MIFno)%DACS(:,DACSno))
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ProcessDACSdata

!=============================================================================
  SUBROUTINE FixLostCarryBits
!=============================================================================

    USE MLSL1Common, ONLY: DACSnum, MaxMIFs
    USE L0_sci_tbls, ONLY: DACS_MAF
    USE Sort_M

    INTEGER :: DACSno, MIFno, Median, BitNo, IntBit1(4), Ndif(4), NumZmatch
    INTEGER :: i, imin(1), mindx(2)
    INTEGER :: Ntot(0:(MaxMIFs-1)), NtotSort(0:(MaxMIFs-1))
    INTEGER :: E                         ! error from median
    INTEGER :: Ebit, EbitVal             ! error to nearest log 2 (bit & value)
    INTEGER, PARAMETER :: E_thold = 48   ! error threshold

    LOGICAL :: Zmatch(4)

    DO DACSno = 1, DACSnum

! Sum D and find median of sum:

       DO MIFno = 0, (MaxMIFs-1)
          Ntot(MIFno) = SUM (DACS_MAF(MIFno)%D(:,DACSno))
       ENDDO

       NtotSort = Ntot
       CALL Sort (NtotSort, 1, MaxMIFs)
       Median = NtotSort(MaxMIFs/2)

! Find error conditions:

       DO MIFno = 0, (MaxMIFs-1)

          IF (Ntot(MIFno) > 0) THEN
             E = Median - Ntot(MIFno)
             IF (E > E_thold) THEN

                Ebit = NINT (LOG (REAL(E)) / LOG (2.0))
                EbitVal = 2**Ebit

! Check lowest bit number set to 1

                Zmatch = .FALSE.   ! no match of trailing Zeros yet
                DO i = 1, 4
                   BitNo = 0
                   IF (DACS_MAF(MIFno)%D(i,DACSno) /= 0) THEN
                      DO
                         IF (BTEST (DACS_MAF(MIFno)%D(i,DACSno), BitNo)) EXIT
                         BitNo = BitNo + 1
                      ENDDO
                   ENDIF
                   IntBit1(i) = 2**BitNo   ! Value of lowest bit set to 1
                   Zmatch(i) = (IntBit1(i) >= EbitVal)   ! match if error
                ENDDO

                NumZmatch = COUNT (Zmatch)

                IF (NumZmatch == 1) THEN   ! Adjust only one counter

                   WHERE (Zmatch)
                      DACS_MAF(MIFno)%D(:,DACSno) = &
                           DACS_MAF(MIFno)%D(:,DACSno) + EbitVal
                   END WHERE

                ELSE IF (NumZmatch > 1) THEN   ! Adjust when multiple counters

                   IF (MIFno == 0) THEN
                      mindx = (/ 1, 1 /)    ! Can't use previous MIF
                   ELSE
                      mindx = (/ -1, 1 /)   ! Use both previous and next MIF
                   ENDIF

                   Ndif = HUGE (Ndif)     ! Init for comparisons
                   WHERE (Zmatch)
                      Ndif = DACS_MAF(MIFno)%D(:,DACSno) - &
                           (DACS_MAF(MIFno+mindx(1))%D(:,DACSno) + &
                           DACS_MAF(MIFno+mindx(2))%D(:,DACSno)) / 2
                   END WHERE
                   imin = MINLOC (Ndif)
                   DACS_MAF(MIFno)%D(imin(1):,DACSno) = &
                        DACS_MAF(MIFno)%D(imin(1):,DACSno) + EbitVal

                ELSE

                   EbitVal = 2**(Ebit-1)   ! Reset value to previous bit
                   Zmatch = .FALSE.
                   DO i = 1, 4
                      Zmatch(i) = (IntBit1(i) >= EbitVal)   ! match if error
                   ENDDO
                   NumZmatch = COUNT (Zmatch)

                   IF (NumZmatch == 2) THEN  ! Only if exactly 2 match
                      WHERE (Zmatch)
                         DACS_MAF(MIFno)%D(:,DACSno) = &
                              DACS_MAF(MIFno)%D(:,DACSno) + EbitVal
                      END WHERE
                   ELSE
                      DACS_MAF(MIFno)%D(:,DACSno) = 0   ! Clear for discard
                   ENDIF
                ENDIF

!!$                write (66, *) 'DACSno, mifno, e, ebitvntot, median: ', &
!!$                     DACSno, MIFno, E, EbitVal, Ntot(MIFno), median
!!$                write (66, *) "D(curr): ", DACS_MAF(MIFno)%D(:,DACSno)
!!$                write (66, *) "D(prev): ", DACS_MAF(MIFno-1)%D(:,DACSno)
!!$                write (66, *) "D(next): ", DACS_MAF(MIFno+1)%D(:,DACSno)
!!$                write (66, "(i6, ':', B20)") (DACS_MAF(MIFno)%D(i,DACSno), &
!!$                     DACS_MAF(MIFno)%D(i,DACSno),i=1,4)
!!$                write (66, *) "IntBit1: ", IntBit1
!!$                write (66, *) "Zmatch: ", Zmatch
!!$                write (66, *) "NumZmatch: ", NumZmatch

             ENDIF
          ENDIF
       ENDDO

    ENDDO


  END SUBROUTINE FixLostCarryBits

!=============================================================================
  SUBROUTINE UnpackDACSdata (C_K, D, Compressed, nchans, R2a)
!=============================================================================

    INTEGER :: C_K(*), D(4), nchans
    LOGICAL :: Compressed
    REAL :: R2a(*)

    REAL :: Ntot3

    IF (Compressed) THEN
       CALL UnpackCompDACS (C_K, R2a)  ! C array type
    ELSE
       Ntot3 = SUM(D) * 3.0
       IF (Ntot3 > 0) CALL UnpackUncompDACS (C_K, Ntot3, nchans, R2a)  ! K array
    ENDIF

  END SUBROUTINE UnpackDACSdata

!=============================================================================
  SUBROUTINE UnpackCompDACS (C, R2a)
!=============================================================================

    USE MLSL1Common, ONLY: DACS_const

    INTEGER :: C(*)
    REAL :: R2a(*)

    INTEGER :: i
    REAL :: L

    L = DACS_const%L    ! convert to REAL

    R2a(1) = 1.0
    DO i = 2, 129
       R2a(i) = (C(i-1) + DACS_const%A(i-1)) / L
    ENDDO

  END SUBROUTINE UnpackCompDACS

!=============================================================================
  SUBROUTINE UnpackUncompDACS (K, Ntot3, nchans, R2a)
!=============================================================================

    INTEGER :: K(*), nchans
    REAL :: R2a(*), Ntot3

    INTEGER :: i

    R2a(1:129) = 0.0
    IF (K(1) == NINT (Ntot3)) RETURN

    R2a(1) = 1.0
    DO i = 2, nchans
       R2a(i) = (K(i) - Ntot3) / (K(1) - Ntot3)
    ENDDO

  END SUBROUTINE UnpackUncompDACS

!=============================================================================
  SUBROUTINE ProcessUnpackedDACS (D, R2a, nchans, TP, DACS_dat)
!=============================================================================

    USE MLSL1Common, ONLY: r8, DACSchans
    USE MathUtils, ONLY: derfi
    
    REAL :: R2a(DACSchans)
    INTEGER :: D(4), TP, nchans
    REAL :: DACS_dat(DACSchans)

    INTEGER :: Dtot, i
    REAL(r8) :: rho(DACSchans), P_thold, N_thold, Z_thold
    REAL(r8) :: A, M, R(256), R_FFT(256)
    REAL(r8), PARAMETER :: sqrt2 = 1.41421356237d0
    REAL(r8), PARAMETER :: C(14) = (/ 0.97523849075787, -0.02380608085090, &
         0.02319499998418, 0.00008254427432, -0.13041555636211, &
         0.07971713086422, 0.00585091297055, -0.06240215483141, &
         0.18410829607653, 0.36609324008142, -0.37590450311119, &
         2.65674351163174, 2.53926887918654, 5.41351505425882 /)

    rho = 0.0
    rho(1) = 1.0

    Dtot = SUM (D)
    IF (Dtot > 0.0) THEN
       P_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * D(1) / Dtot)
       N_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * D(4) / Dtot)
       Z_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * (D(1) + D(2))/ Dtot)
    ELSE
       P_thold = 0.0
       N_thold = 0.0
       Z_thold = 0.0
    ENDIF
    A = P_thold - N_thold
    M = (P_thold + N_thold) / 2.0 - 0.9

    DO i = 2, nchans
       rho(i) = C(1) * R2a(i) + C(2) * R2a(i)**3 + C(3) * R2a(i)**5 + &
            C(4) * R2a(i)**7 + C(5) * SIN(R2a(i) * C(12)) * M + &
            C(6) * SIN(R2a(i) * C(13)) * M**2 + C(7) * SIN(R2a(i)*C(14)) * M + &
            C(8) * A**2 + C(9) * A**2 * R2a(i) + C(10) * Z_thold**2*R2a(i) + &
            C(11) * Z_thold * A
    ENDDO

! Apodize here?

! Prepare for FFT

  R(1:129) = rho
  R(130:256) = rho(128:2:-1)

! Do the FFT

  CALL rfftw_f77_one (plan129, R, R_FFT)

  DACS_dat = R_FFT(1:129) * TP

  END SUBROUTINE ProcessUnpackedDACS

END MODULE DACsUtils

! $Log$
! Revision 2.4  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.3  2003/03/25 19:53:54  perun
! Test D before bit test
!
! Revision 2.2  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.1  2002/03/29 20:20:16  perun
! Version 1.0 commit
!
!
