! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE DACsUtils
!=============================================================================

  USE MLSL1Utils, ONLY: BigEndianStr
  USE L1BData, ONLY: L1BData_T, ReadL1BData, DeallocateL1BData

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ExtractDACSdata, UncompressDACSdata, ProcessDACSdata, InitDACS_FFT
  PUBLIC :: FinalizeDACSdata

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
    INTEGER, DIMENSION(:), INTENT(OUT) :: C
    INTEGER, INTENT(OUT) :: D(4), TP, DIO, Zlag, LO

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

             ENDIF
          ENDIF
       ENDDO

    ENDDO


  END SUBROUTINE FixLostCarryBits

!=============================================================================
  SUBROUTINE UnpackDACSdata (C_K, D, Compressed, nchans, R2a)
!=============================================================================

    INTEGER, DIMENSION(:) :: C_K
    INTEGER :: D(4), nchans
    LOGICAL :: Compressed
    REAL, DIMENSION(:) :: R2a

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

    INTEGER, DIMENSION(:) :: C
    REAL, DIMENSION(:) :: R2a

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

    INTEGER, DIMENSION(:) :: K
    INTEGER :: nchans
    REAL, DIMENSION(:) :: R2a
    REAL :: Ntot3

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

    IF (D(1) == 0 .AND. D(4) == 0) THEN
       DACS_dat = 0.0
       RETURN
    ENDIF

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

!=============================================================================
  SUBROUTINE FinalizeDACSdata
!=============================================================================

    USE MLSL1Common, ONLY: L1BFileInfo, DACSchans, DACSnum
    USE MLSCommon, ONLY: rm
    USE MLSL1Config, ONLY: MIFsGHz, L1Config
    USE MLS_DataProducts, ONLY: DataProducts_T, Deallocate_DataProducts
    USE MLSAuxData, ONLY: Build_MLSAuxData
    USE MLSHDF5, ONLY: IsHDF5DSPresent, MakeHDF5Attribute
    USE HDF5, ONLY: h5gopen_f
    USE MatrixModule_0, ONLY: MatrixInversion
    USE units, ONLY: Pi

    INTEGER :: bno, ch, grp_id, i, Flag, maf, mif, noMAFs, sd_id, status
    INTEGER, PARAMETER :: rch(2) = (/ 10, 100 /)  ! range of channels
    INTEGER, PARAMETER :: Nch = (rch(2) - rch(1) + 1)  ! number of "good" chans
    INTEGER :: Ycount(nch)
    REAL :: Y(nch)  ! will contain average over + prec and within alt range
    REAL, POINTER, DIMENSION(:,:) :: alt, sza
    REAL, POINTER, DIMENSION(:,:,:) :: rad, rad_prec
    REAL, PARAMETER :: MinAlt = 78000, MaxAlt = 100000 ! Altitude ranges (m)
 
! High altitude model parameters:

    REAL, PARAMETER :: cen(24:25) = (/ 48.384, 51.44 /)
    REAL, PARAMETER :: dopp(24:25) = (/ 1.667, 2.2 /)
    REAL, PARAMETER :: spur = (5.0e6 / 7*128 / 1.25e7)
    REAL, PARAMETER :: f(Nch) = (/ ((i-1), i=rch(1),rch(2)) /)
    REAL(rm) :: X(4,Nch), Xinv(4,4), Avec(4,24:25), B(4,Nch)
    REAL :: apod(DACSchans,24:25), spurmag(2,24:25)
    REAL :: NoiseInflationFactor(24:25)
    REAL, PARAMETER :: MaxNoiseInflationFactor = 3.0
    REAL, PARAMETER :: MinSZA(24:25) = (/ 100.0, 0.0 /)

    TYPE (L1BData_T) :: L1BData
    TYPE (DataProducts_T) :: DACsDS  

    CHARACTER(len=26), PARAMETER :: DACS_Name(24:25) = (/ &
         "R3:240.B24D:O3.S0.DACS-3  ", "R3:240.B25D:CO.S1.DACS-1  " /)

!=============================================================================

    IF (.NOT. L1Config%Calib%CalibDACS) RETURN   ! Nothing to do

    PRINT *, 'Finalizing DACS data...'

! Initialize attributes:

    apod = 0.0
    Avec = 0.0
    spurmag = 0.0
    NoiseInflationFactor = 0.0

! Get altitude and SZA (GHz):

    sd_id = L1BFileInfo%OAid

    CALL ReadL1BData (sd_id, '/GHz/GeodAlt', L1BData, noMAFs, Flag, &
         NeverFail=.TRUE., HDFversion=5)
    ALLOCATE (alt(MIFsGHz,noMAFs))
    alt = L1BData%DpField(1,:,:)
    CALL DeallocateL1BData (L1BData)

    CALL ReadL1BData (sd_id, '/GHz/SolarZenith', L1BData, noMAFs, Flag, &
         NeverFail=.TRUE., HDFversion=5)
    ALLOCATE (sza(MIFsGHz,noMAFs))
    sza = L1BData%DpField(1,:,:)
    CALL DeallocateL1BData (L1BData)

! Set up for HDF output:

    CALL Deallocate_DataProducts (DACsDS)
    ALLOCATE (DACsDS%Dimensions(3))
    DACsDS%data_type = 'real'
    DACsDS%Dimensions(1) = 'chanDACS'
    DACsDS%Dimensions(2) = 'GHz.MIF             '
    DACsDS%Dimensions(3) = 'MAF                 '

! get Band 24/25 data:

    sd_id = L1BFileInfo%RADDid

    DO bno = 25, 24, -1

       IF (.NOT. IsHDF5DSPresent (sd_id, TRIM(DACS_Name(bno)))) CYCLE

       DEALLOCATE (rad, stat=status)
       DEALLOCATE (rad_prec, stat=status)
       CALL ReadL1BData (sd_id, DACS_Name(bno), L1BData, noMAFs, Flag, &
            NeverFail=.TRUE., HDFversion=5)
       ALLOCATE (rad(DACSchans,MIFsGHz,noMAFs))
       rad = L1BData%DpField(:,:,:)
       CALL ReadL1BData (sd_id, TRIM(DACS_Name(bno))//' precision', L1BData, &
            noMAFs, Flag, NeverFail=.TRUE., HDFversion=5)
       ALLOCATE (rad_prec(DACSchans,MIFsGHz,noMAFs))
       rad_prec = L1BData%DpField(:,:,:)
       CALL DeallocateL1BData (L1BData)

! Overwrite channels 110:129 with average of channels 100:110:

       DO maf = 1, noMAFs
          DO mif = 1, MIFsGHz
             rad(110:129,mif,maf) = SUM(rad(100:110,mif,maf)) / 11
          ENDDO
       ENDDO
       rad_prec(110:129,:,:) = -1 * ABS(rad_prec(110:129,:,:))  ! Negate precs

! Mean over good, high altitude data:

       Y = 0.0
       Ycount = 0
       DO maf = 1, noMAFs
          DO mif = 1, MIFsGHz
             IF (alt(mif,maf) >= MinAlt .AND. alt(mif,maf) <= MaxAlt .AND. &
                  sza(mif,maf) > MinSZA(bno)) THEN
                DO ch = rch(1), rch(2)  ! channel range
                   IF (rad_prec(ch,mif,maf) > 0.0) THEN
                      i = ch - rch(1) + 1
                      Y(i) = Y(i) + rad(ch,mif,maf)
                      yCount(i) = yCount(i) + 1
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO

       DO ch = 1, Nch
          IF (Ycount(ch) > 0) Y(ch) = Y(ch) / Ycount(ch)
       ENDDO

! High altitude model radiances:

       X(1,:) = 2*EXP(-((f-cen(bno))/dopp(bno))**2/2)
       X(2,:) = EXP(-((f-cen(bno)-spur)/dopp(bno))**2/2) + &
            EXP(-((f-cen(bno)+spur)/dopp(bno))**2/2)
       X(3,:) = EXP(-((f-cen(bno)-2*spur)/dopp(bno))**2/2) + &
            EXP(-((f-cen(bno)+2*spur)/dopp(bno))**2/2)
       X(4,:) = 1.0
       Xinv = MATMUL (X, TRANSPOSE(X))

       CALL MatrixInversion (Xinv)

       B = MATMUL (Xinv, X)
       Avec(:,bno) = MATMUL (B, Y)

       spurmag(1,bno) = Avec(2,bno) / SUM (Avec(1:3,bno))
       spurmag(2,bno) = Avec(3,bno) / SUM (Avec(1:3,bno))

       DO i = 1, DACSchans
          apod(i,bno) = 1.0 / (1 + spurmag(1,bno)*(COS(5e6/7*Pi/1.25e7 * &
               (i-1)) - 1) + spurmag(2,bno)*(COS(2*5e6/7*Pi/1.25e7 * (i-1))-1))
       ENDDO

       NoiseInflationFactor(bno) = SQRT (SUM(apod(:,bno)**2) / 129)

       IF (NoiseInflationFactor(bno) < MaxNoiseInflationFactor) THEN

! Deconvolve:

          CALL DeconvolveRads (rad, apod(:,bno))

! Scale precisions:

          rad_prec = rad_prec * NoiseInflationFactor(bno)

! Output Band rads:

          DACsDS%name = DACS_Name(bno)
          DO maf = 1, noMAFS
             CALL Build_MLSAuxData (sd_id, DACsDS, rad(:,:,maf), &
                  lastIndex=maf, disable_attrib=.TRUE.)
          ENDDO

       ELSE
          rad_prec = -1 * ABS(rad_prec)  ! Negate precs
       ENDIF

! Output Band precisions:

       DACsDS%name = TRIM(DACS_Name(bno))//' precision'
       DO maf = 1, noMAFS
          CALL Build_MLSAuxData (sd_id, DACsDS, rad_prec(:,:,maf), &
               lastIndex=maf, disable_attrib=.TRUE.)
       ENDDO

    ENDDO

! Output attribute vectors:

    CALL h5gopen_f (sd_id, '/', grp_id, status)
    CALL MakeHDF5Attribute (grp_id, 'apodB24', apod(:,24), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'apodB25', apod(:,25), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'AvecB24', Avec(:,24), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'AvecB25', Avec(:,25), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'spurmagB24', spurmag(:,24), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'spurmagB25', spurmag(:,25), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'NoiseInflationFactorB24', &
         NoiseInflationFactor(24:24), .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'NoiseInflationFactorB25', &
         NoiseInflationFactor(25:25), .TRUE.)

  END SUBROUTINE FinalizeDACSdata

!=============================================================================
  SUBROUTINE DeconvolveRads (rad, apod)
!=============================================================================

    USE MLSCommon, ONLY: r8
    USE DFFT_m, ONLY: DRFT1

    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: rad
    REAL, DIMENSION(:), INTENT(IN) :: apod

    REAL(r8) :: dacs_dat(256), S(256)
    INTEGER :: ms, maf, MAFs, mif, MIFs

    MIFs = SIZE (rad(1,:,1))
    MAFs = SIZE (rad(1,1,:))

    ms = 0
    DO mif = 1, MIFs
       DO maf = 1, MAFs
          dacs_dat(1:129) = rad(:,mif,maf)
          dacs_dat(130:256) = dacs_dat(128:2:-1)

! Analyze:

          CALL DRFT1 (dacs_dat, 'A', 8, ms, S)

! Apodize:

          dacs_dat(1) = dacs_dat(1) * apod(1)
          dacs_dat(2) = dacs_dat(2) * apod(129)
          dacs_dat(3:255:2) = dacs_dat(3:255:2) * apod(2:128)
          dacs_dat(4:256:2) = dacs_dat(4:256:2) * apod(2:128)

! Synthesize:

          CALL DRFT1 (dacs_dat, 'S', 8, ms, S)
          rad(:,mif,maf) = dacs_dat(1:129)

       ENDDO
    ENDDO

  END SUBROUTINE DeconvolveRads

END MODULE DACsUtils

! $Log$
! Revision 2.8  2005/01/25 15:19:20  perun
! Add TRIM to remove possible extra blank in precision dataset
!
! Revision 2.7  2004/12/01 17:08:56  perun
! Add routines to deconvolve and remove spurs
!
! Revision 2.6  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.5  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
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
