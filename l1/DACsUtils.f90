! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE DACsUtils
!=============================================================================

  USE MLSL1Utils, ONLY: BigEndianStr

  IMPLICIT NONE

  INTEGER :: plan129

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE InitDACS_FFT

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

  SUBROUTINE ExtractDACSdata (rawdata, K, D, TP, DIO, LO, nchans)

    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: rawdata
    INTEGER, DIMENSION(:), INTENT(OUT) :: K
    INTEGER, INTENT(OUT) :: D(4), TP, DIO, LO, nchans

    INTEGER :: i, offset

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

  SUBROUTINE UncompressDACSdata (rawdata)

    CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: rawdata

    INTEGER :: C(128), D(4), TP, DIO, Zlag, LO
    INTEGER :: i, i12, i8, rindx

!! Masks for extending sign bits:

    INTEGER, PARAMETER :: bitmask7 = z'80'
    INTEGER, PARAMETER :: signmask7 = z'FFFFFF00'
    INTEGER, PARAMETER :: bitmask11 = z'800'
    INTEGER, PARAMETER :: signmask11 = z'FFFFF000'

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
       i12 = IAND (i12, z'FFF')        ! mask off upper 4 bits
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

  SUBROUTINE ProcessDACSdata (D, DACS_K, nchans, TP, DACS_dat)

    USE MLSL1Common, ONLY: r8, DACSchans
    USE MathUtils, ONLY: derfi
    
    INTEGER :: DACS_K(DACSchans)
    INTEGER :: D(4), TP, nchans
    REAL :: DACS_dat(DACSchans)

    INTEGER :: Dtot, DIO, LO, i
    REAL(r8) :: lag0, offset, rho_hat, rho(DACSchans), P_thold, N_thold, Z_thold
    REAL(r8) :: A, M, R(256), R_FFT(256)
    REAL(r8), PARAMETER :: sqrt2 = 1.41421356237d0
    REAL(r8), PARAMETER :: C(14) = (/ 0.97523849075787, -0.02380608085090, &
         0.02319499998418, 0.00008254427432, -0.13041555636211, &
         0.07971713086422, 0.00585091297055, -0.06240215483141, &
         0.18410829607653, 0.36609324008142, -0.37590450311119, &
         2.65674351163174, 2.53926887918654, 5.41351505425882 /)

    offset = 3.0 * SUM (D)
    rho = 0.0
    rho(1) = 1.0
    lag0 = DACS_K(1) - offset

    Dtot = SUM (D)
    IF (Dtot > 0.0) THEN
       P_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * D(1) / Dtot)
       N_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * D(4) / Dtot)
       Z_thold = sqrt2 * derfi (1.0d0 - 2.0d0 * (D(1) + D(2))/ Dtot)
!PRINT *, "P/N/Z_thold: ", REAL(P_thold), REAL(N_thold), REAL(Z_thold)
    ELSE
       P_thold = 0.0
       N_thold = 0.0
       Z_thold = 0.0
    ENDIF
    A = P_thold - N_thold
    M = (P_thold + N_thold) / 2.0 - 0.9

    DO i = 2, nchans
       IF (lag0 /= 0.0) THEN
          rho_hat = (DACS_K(i) - offset) / lag0
       ELSE
          rho_hat = 0.0
       ENDIF
       rho(i) = C(1) * rho_hat + C(2) * rho_hat**3 + C(3) * rho_hat**5 + &
            C(4) * rho_hat**7 + C(5) * SIN(rho_hat * C(12)) * M + &
            C(6) * SIN(rho_hat * C(13)) * M**2 &
            + C(7) * SIN(rho_hat * C(14)) * M + C(8) * A**2 + &
            C(9) * A**2 * rho_hat + C(10) * Z_thold**2 * rho_hat + &
            C(11) * Z_thold * A
    ENDDO

! Apodize here?


!  rho = rho * TP
! print *, 'rho: ', rho

! Prepare for FFT

  R(1:129) = rho
  R(130:256) = rho(128:2:-1)

! Do the FFT

  CALL rfftw_f77_one (plan129, R, R_FFT)

  DACS_dat = R_FFT(1:129) * TP
!!$    if (d(1) == 48303 .and. d(2) == 81510) then
!!$       print *, DACS_dat
!!$    endif

  END SUBROUTINE ProcessDACSdata

END MODULE DACsUtils

! $Log$
! Revision 2.1  2002/03/29 20:20:16  perun
! Version 1.0 commit
!
!
