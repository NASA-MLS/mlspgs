PROGRAM conv_uars

USE rad_file_contents
USE swapendian
USE MLSFiles, ONLY: MLS_openFile, MLS_closeFile
USE MLSAuxData, ONLY: CreateGroup_MLSAuxData
USE MLSHDF5, ONLY: MLS_h5open, MLS_h5close
USE Output_UARS_L1B
USE SDPToolkit
USE oa_file_contents

IMPLICIT NONE

CHARACTER (len=80) :: arg, infile='uars.dat', Rad_File='RAD.h5', &
     OA_File='OA.h5', out_dir='/data/'
INTEGER, PARAMETER :: in_unit=101, inlen=14336/4, HDFversion=5
INTEGER :: i, j, k, mafno=0, nummmaf, recno, rad_sd_id, error, oa_sd_id, stat
INTEGER :: year, doy, month, day, ms, hrs, mins, secs, yrdoy, noargs
INTEGER, DIMENSION(32) :: mifno = (/ (i, i=1,32) /)

REAL :: allrad(15*6*32), rad(15*6,32), rad6(15,6,32), &
     allprec(15*6*32), prec(15*6,32), prec6(15,6,32)
EQUIVALENCE (allrad, rad), (allprec, prec), (rad, rad6), (prec, prec6)
CHARACTER (LEN=*), PARAMETER :: &
     DateFmt = "(I4, '-', I3.3, 'T', I2.2, ':', I2.2, ':', I2.2)", &
     data_dir = '/data/umls/l1/'
CHARACTER (len=17) :: utc_time
CHARACTER (len=8) :: yrdoy_str = 'yyyydnnn'
REAL*8 :: tai_time, gird_time
REAL*8, PARAMETER :: TAItoGIRD = 1104537627.0d00

INTEGER*2, EXTERNAL :: SwapShort  ! should put into SwapEndian later!!!

! Source code:

! Open UARS input file:

noargs = COMMAND_ARGUMENT_COUNT()
IF (noargs == 0) THEN
   PRINT *, 'No input file provided.'
   PRINT *, 'Using defaults...'
ELSE
   i = 1

   DO WHILE (i <= noargs)

      CALL GET_COMMAND_ARGUMENT (i, arg)

      SELECT CASE (arg)
      CASE default
         !infile = data_dir // TRIM(arg)     ! input filename (with dir)
         infile = TRIM(arg)                 ! input filename
         IF (noargs /= 3) THEN
            PRINT *, 'Missing output directory argument!'
            CALL EXIT
         ENDIF
      CASE ('-o', '--outdir')
         i = i + 1
         CALL GET_COMMAND_ARGUMENT (i, out_dir)  ! next argument for out dir
      CASE ('-h', '--help')
         PRINT '(a)', 'usage: conv_uars [-h] input_file -o output_dir' 
         CALL EXIT
      END SELECT
      i = i + 1   ! Next argument
   ENDDO
ENDIF

OPEN (unit=in_unit, file=infile, status="OLD", form="UNFORMATTED", &
     access="DIRECT", recl=inlen, iostat=error)

IF (error /= 0) THEN
   PRINT *, 'No such input file: '//infile
   CALL EXIT
ENDIF

! read and convert header part:

DO recno = 1, 3

   READ (unit=in_unit, rec=recno) lvl1_hdr

   ! convert to little endian as needed:

   lvl1_hdr%recordno = SwapBig (lvl1_hdr%recordno)
   DO i = 1, 2
      lvl1_hdr%write_time(i) = SwapBig (lvl1_hdr%write_time(i))
      lvl1_hdr%start_time(i) = SwapBig (lvl1_hdr%start_time(i))
      lvl1_hdr%stop_time(i) = SwapBig (lvl1_hdr%stop_time(i))
   ENDDO
   lvl1_hdr%nummmaf = SwapBig (lvl1_hdr%nummmaf)
   lvl1_hdr%uars_day = SwapBig (lvl1_hdr%uars_day)
   lvl1_hdr%mls_status_day = SwapBig (lvl1_hdr%mls_status_day)

ENDDO

nummmaf = lvl1_hdr%nummmaf

IF (nummmaf == 0) THEN
   PRINT *, 'No data in file: '//infile
   CALL EXIT
ENDIF

! get the yrdoy of file:

recno = nummmaf / 2
READ (unit=in_unit, rec=recno) limb_hdr
yrdoy = SwapBig (limb_hdr%mmaf_time(1))

IF (noargs > 0) THEN   ! use passed input for setting outputs
   year = yrdoy / 1000 + 1900
   doy = MOD (yrdoy, 1000)
   WRITE (yrdoy_str, '(i4.4, "d", i3.3)') year, doy
   PRINT *, 'yrdoy, year, doy, str: ', yrdoy, year, doy, yrdoy_str

   Rad_File = TRIM(out_dir) // '/MLS-UARS_L1BRAD_' // yrdoy_str // '.h5'
   OA_File = TRIM(out_dir) // '/MLS-UARS_L1BOA_' // yrdoy_str // '.h5'
ENDIF

PRINT *, 'infile: ', infile
PRINT *, 'rad file: ', Rad_File
PRINT *, 'OA file: ', OA_File

! Open output RAD and OA HDF files:

error = 0
CALL MLS_h5open (error)
 
CALL MLS_openFile (Rad_File, 'create', rad_sd_id, HDFversion)

CALL MLS_openFile (OA_File, 'create', oa_sd_id, HDFversion)
CALL CreateGroup_MLSAuxData (oa_sd_id, 'GHz')
CALL CreateGroup_MLSAuxData (oa_sd_id, 'sc')

! save attributes for both HDF files:

CALL set_attributes (rad_sd_id, yrdoy, lvl1_hdr%start_time, lvl1_hdr%stop_time)
CALL set_attributes (oa_sd_id, yrdoy, lvl1_hdr%start_time, lvl1_hdr%stop_time)

! print stuff!

PRINT *, 'uars day', lvl1_hdr%uars_day
PRINT *, 'start time: ', lvl1_hdr%start_time
PRINT *, 'stop time: ', lvl1_hdr%stop_time
PRINT *, 'yrdoy: ', yrdoy

PRINT *, 'nummmaf: ', nummmaf
! read and convert data:

DO recno = 4, nummmaf+3

   READ (unit=in_unit, rec=recno) limb_hdr, limb_stat, limb_rad, limb_oa
   mafno = mafno + 1

   ! convert to little endian as needed:

   limb_hdr%recordno = SwapBig (limb_hdr%recordno)
   limb_hdr%mmafno = SwapBig (limb_hdr%mmafno)
   DO i = 1, 2
      limb_hdr%mmaf_time(i) = SwapBig (limb_hdr%mmaf_time(i))
   ENDDO
   DO i = 1, 16
      limb_stat%prd_temps(i) = SwapBig (limb_stat%prd_temps(i))
   ENDDO
   limb_stat%maneuver_stat = SwapBig (limb_stat%maneuver_stat)
   limb_stat%mls_status = SwapBig (limb_stat%mls_status)
   DO i = 1, 90
       limb_stat%wall_mmaf(i) = SwapShort (limb_stat%wall_mmaf(i))
   ENDDO

   ! convert radiances:

   DO k = 1, 32           ! MIF number
      DO j = 1, 2         ! radiance/error
         DO i = 1, 90     ! radiometer and channel
            limb_rad%rad_l1(i,j,k) = SwapShort (limb_rad%rad_l1(i,j,k))
         ENDDO
      ENDDO
   ENDDO

   rad = limb_rad%rad_l1(:,1,:)
   prec = limb_rad%rad_l1(:,2,:)
   WHERE (allrad .LT. 0.0 .AND. allrad .NE. -1000.0) allrad = allrad + 65536.0
   WHERE (allprec .LT. 0.0 .AND. allprec .NE. -1000.0) &
        allprec = allprec + 65536.0
   WHERE (allrad .NE. -1000.0) allrad = allrad * 0.01   ! scale radiances

   WHERE (allprec .EQ. -1000.0) allprec = -1.0
   WHERE (allprec .NE. -1.0) allprec = allprec * 0.01   ! scale precisions

   year = limb_hdr%mmaf_time(1) / 1000 + 1900
   doy = MOD (limb_hdr%mmaf_time(1), 1000)

   ! convert millisecs to HMS:

   CALL ms_to_hms (limb_hdr%mmaf_time(2), hrs, mins, secs)
   WRITE (utc_time, fmt=DateFmt) year, doy, hrs, mins, secs

   ! get and save the GIRD time:

   stat = PGS_TD_UTCtoTAI (utc_time, tai_time)
   gird_time = tai_time + TAItoGIRD

   ! write radiances etc. to file:

   CALL OutputL1B_rad (rad_sd_id, mafno, limb_hdr%mmafno, rad6, prec6, &
                       gird_time)

   ! swap all OA bytes:

   CALL swap_oa_rec

   ! convert to EMLS OA:

   CALL uars_to_emls_oa

   ! write to OA file:

   CALL OutputL1B_OA (oa_sd_id, mafno, limb_hdr%mmafno)

ENDDO

PRINT *, 'mafs: ', mafno 

CLOSE (in_unit)

CALL MLS_closeFile (rad_sd_id, HDFversion=HDFversion)
CALL MLS_closeFile (oa_sd_id, HDFversion=HDFversion)
CALL MLS_h5close (error)

END PROGRAM conv_uars
