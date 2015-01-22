! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Conv_UARS

  use Constants, only: Deg2Rad
use Constants, only: Rad2Deg
  use HDF5, only:  H5GClose_f, H5GOpen_f
  use ISO_Fortran_Env, only: Output_Unit
  use MLSAuxData, ONLY: CreateGroup_MLSAuxData
  use MLSCommon, ONLY: FileNameLen
  use MLSFiles, ONLY: MLS_openFile, MLS_closeFile
  use MLSHDF5, ONLY: MakeHDF5Attribute, MLS_h5open, MLS_h5close
  use MLSKinds, only: R8
  use Output_UARS_L1B, only: OutputL1B_OA, OutputL1B_Rad
  use Rad_File_Contents, only: Lvl1_Hdr_t, Limb_Hdr_t, Limb_Stat_t, Limb_Rad_t, &
                             & Limb_OA_t
  use Set_Attributes_m, only: Set_Attributes
  use SwapEndian, only: SwapBig
  use Swap_OA_Rec_m, only: Swap_OA_Rec
  use UARS_to_EMLS_OA_m, only: EarthEllipseTag, Frame_TAI, NoMIFs

  implicit NONE

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Is the input file name stored as an attribute in the output file?
  logical, parameter :: Input_File_Attrib = .true.

  integer, parameter :: In_Unit=101, InLen=14336/4, HDFversion=5
  character (len=FileNameLen) :: Arg
  character (len=FileNameLen) :: InFile='uars.dat'
  character (len=FileNameLen) :: OA_File='OA.h5'
  character (len=FileNameLen) :: Out_Dir='/data/'
  character (len=FileNameLen) :: Rad_File='RAD.h5'
  integer :: DOY
  integer :: Error
  integer :: Grp_ID  ! of root directory of output file
  integer :: I
  integer :: MAF1    ! First MAF, default 1
  integer :: MAFn    ! Last MAF, default from header
  integer :: MAFno
  integer :: n_days = 0
  integer :: NoArgs
  integer :: NummMAF
  integer :: OA_SD_id
  integer :: Rad_SD_id
  integer :: Recno
  integer :: Stat
  integer :: Year
  integer :: YrDOY
  integer :: W ! zero or one, to select Limb_Hdr etc.

  ! Inputs
  type(lvl1_hdr_t) :: Lvl1_Hdr
  type(limb_hdr_t) :: Limb_Hdr(0:1)
  type(limb_stat_t) :: Limb_Stat(0:1)
  type(limb_rad_t) :: Limb_Rad(0:1)
  type(limb_oa_t) :: Limb_OA(0:1)

  character (len=25) :: AsciiUTC(0:1)
  character (len=8) :: n_str
  character (len=17) :: UTC_Time
  character (len=8) :: YrDoy_str = 'yyyydnnn'
  real(r8) :: GIRD_Time
  real(r8) :: SecTAI93
  real(r8) :: Sat_Vel(3,0:1), vN(0:1)  ! SC Velocity in ECR, its length, km/s
  real(r8) :: Vel_ECI(3,0:1), vNI(0:1) ! SC Velocity in ECR, its length, km/s
  real(r8) :: PosECR(3,0:1)            ! SC Position, in ECR, meters
  real(r8) :: dVV(3)                   ! delta V in ECR
  real(r8) :: dV, dVN                  ! delta V angle, |delta V|, ECR, per MAF
  real(r8) :: dVI, dVIN                ! delta V angle, |delta V|, ECI, per MAF
  real(r8) :: ECI(6,1), ECR(6,1)       ! Position, velocity, for conversions
  real(r8) :: TAI_Time
  real(r8), parameter :: TAItoGIRD = 1104537627.0d00

  ! Source code:
  ! Don't wrap stdout at 80 chars
  open ( unit=output_unit, recl=4096 )

  ! Open UARS input file:

  MAF1 = 1
  MAFn = -1  ! Signal to get it from the header

  noargs = COMMAND_ARGUMENT_COUNT()
  if (noargs == 0) then
    print *, 'No input file provided.'
    print *, 'Using defaults...'
  else
    i = 1

    do while (i <= noargs)

      call GET_COMMAND_ARGUMENT (i, arg)

      select case (arg)
      case default
        infile = TRIM(arg)                 ! input filename
        if (noargs < 3) then
          print *, 'Missing output directory argument? Too few args!'
          go to 9
        end if
      case ('-b', '--backdate')
        i = i + 1
        call GET_COMMAND_ARGUMENT (i, n_str)  ! how many days to backdate times in file
        read( n_str, * ) n_days
        ! Now we must take account for how mmaf_time actually encodes dates
        ! It uses a convention that the year contains 1000 days
        ! if the year y = 1900 + p
        ! and the doy   = d
        ! then mmaf_time = 1000*y + d (yes, I know this is ugh--
        ! "It was like that when I got here!")
        ! So we're assuming you either want to backdate by one day or so 
        ! or else by a number of years
        ! E.g., to backdate 10 years, set n_days to 3650 (we ignore leap years)
        if ( n_days > 364 ) n_days = 1000 * (n_days/365)
      case ( '-M', '--MAF' )
        i = i + 1
        call GET_COMMAND_ARGUMENT (i, n_str) ! First MAF
        read ( n_str, * ) MAF1
        i = i + 1
        call GET_COMMAND_ARGUMENT (i, n_str) ! Last MAF
        read ( n_str, * ) MAFn
      case ('-o', '--outdir')
        i = i + 1
        call GET_COMMAND_ARGUMENT (i, out_dir)  ! next argument for out dir
      case ('-h', '--help')
        call GET_COMMAND_ARGUMENT ( 0, arg )
        print '(a)', 'usage: ' // trim(arg) // ' [options] input_file' 
        print '(a)', '      options' 
        print '(a)', '  -b n_days         backdate converted files in time by n_days' 
        print '(a)', '  -M MAF1 MAFn      MAF range to process'
        print '(a)', '  -o output_dir     write converted files to output_dir' 
        print '(a)', '  -h                print help' 
        go to 9
      end select
      i = i + 1   ! Next argument
    end do
  end if

  ! The length computed by
  ! inquire ( iolength=i ) limb_hdr(1), limb_stat(1), limb_rad(1), limb_oa(1)
  ! is 90 less than inlen.  What do the other 360 bytes contain?  It looks
  ! like zeroes, but is it just padding or do some records contain
  ! something useful there?

  open ( unit=in_unit, file=infile, status="OLD", form="UNFORMATTED", &
       & access="DIRECT", recl=inlen, iostat=error )

  if (error /= 0) then
    print *, 'No such input file: '//infile
    go to 9
  end if

  ! read and convert header part:

  do recno = 1, 3

    read (unit=in_unit, rec=recno) lvl1_hdr

    ! convert to little endian as needed:

    lvl1_hdr%recordno = SwapBig (lvl1_hdr%recordno)
    do i = 1, 2
      lvl1_hdr%write_time(i) = SwapBig (lvl1_hdr%write_time(i))
      lvl1_hdr%start_time(i) = SwapBig (lvl1_hdr%start_time(i))
      lvl1_hdr%stop_time(i) = SwapBig (lvl1_hdr%stop_time(i))
    end do
    lvl1_hdr%nummmaf = SwapBig (lvl1_hdr%nummmaf)
    lvl1_hdr%uars_day = SwapBig (lvl1_hdr%uars_day)
    lvl1_hdr%mls_status_day = SwapBig (lvl1_hdr%mls_status_day)

  end do

  nummmaf = lvl1_hdr%nummmaf
  if ( MAFn < 0 ) MAFn = nummmaf
  if ( MAF1 < 1 .or. MAFn > nummmaf .or. MAF1 > MAFn) then
    print '(3(a,i0))', 'MAF numbers out of range.  MAF1 = ', MAF1, &
      & ' MAFn = ', MAFn, ' Valid range is 1 .. ', nummmaf
    stop
  end if

  if (nummmaf == 0) then
    print *, 'No data in file: '//infile
    go to 9
  end if

  ! get the yrdoy of file:

  recno = nummmaf / 2
  read ( unit=in_unit, rec=recno ) limb_hdr(0)
  yrdoy = SwapBig (limb_hdr(0)%mmaf_time(1)) - n_days ! This lets us backdate the files; e.g., for sids

  if ( noargs > 0 ) then   ! use passed input for setting outputs
    year = yrdoy / 1000 + 1900
    doy = mod (yrdoy, 1000)
    write (yrdoy_str, '(i4.4, "d", i3.3)') year, doy
    print *, 'yrdoy, year, doy, str: ', yrdoy, year, doy, yrdoy_str

    Rad_File = TRIM(out_dir) // '/MLS-UARS_L1BRAD_' // yrdoy_str // '.h5'
    OA_File = TRIM(out_dir) // '/MLS-UARS_L1BOA_' // yrdoy_str // '.h5'
  end if

  write( 6, '(a)' ) 'infile: ' // trim(infile)
  write( 6, '(a)' ) 'rad file: ' // trim(Rad_File)
  write( 6, '(a)' ) 'OA file: ' // trim(OA_File)
  print *, 'n days backdated: ', n_days

  ! Open output RAD and OA HDF files:

  error = 0
  call MLS_h5open ( error )

  call MLS_openFile ( Rad_File, 'create', rad_sd_id, HDFversion )

  call MLS_openFile ( OA_File, 'create', oa_sd_id, HDFversion )
  call CreateGroup_MLSAuxData ( oa_sd_id, 'GHz' )
  call CreateGroup_MLSAuxData ( oa_sd_id, 'sc' )

  ! save attributes for both HDF files:

  call set_attributes ( rad_sd_id, yrdoy, lvl1_hdr%start_time, lvl1_hdr%stop_time )
  call set_attributes ( oa_sd_id, yrdoy, lvl1_hdr%start_time, lvl1_hdr%stop_time )

  if ( input_file_attrib ) then
    call h5gopen_f ( oa_sd_id, '/', grp_id, stat )
    call MakeHDF5Attribute(grp_id, 'InputFileName', trim(infile), .true.)
  end if

  ! print stuff

  print *, 'UARS day', lvl1_hdr%uars_day
  print *, 'Start time: ', lvl1_hdr%start_time
  print *, 'Stop time:  ', lvl1_hdr%stop_time
  print *, 'yrdoy: ', yrdoy

  print *, 'nummmaf: ', nummmaf
  ! read and convert data:

  MAFno = MAF1 - 1

  w = 0
  read ( unit=in_unit, rec=MAF1+3 ) &
    & limb_hdr(w), limb_stat(w), limb_rad(w), limb_oa(w)
  ! swap all OA bytes:
  call swap_OA_rec ( limb_oa(w) )

  vel_ECI(:,w) = limb_oa(w)%sat_vel
  vNI(w) = norm2(vel_ECI(:,w))

  call sat_vel_in_ecr ( limb_oa(w), asciiUTC(w), sat_vel(:,w), posecr(:,w) )

  vN(w) = norm2(sat_vel(:,w))
  
  dV = 0.0
  dVI = 0.0
  dVN = 0.0
  dVIN = 0.0
  dVV = 0.0
write ( 42,'(a/2x,3f9.2,1p,4g15.6)' ) &
& '     VelECI   VelECR      vNI ----- dV ----- ----- dVN ---- ----- dVI ---- ---- dVIN ----', &
& vN(w), vNI(w), norm2(posecr(:,w))/1000, dV*rad2deg, dVN*1.0e6, dVI*rad2deg, dVIN*1.0e6

  do recno = MAF1+3, MAFn+3

    ! Get the record for the next MAF, if there is one.
    ! Calculate the change in velocity vector angle.
    ! This is used to rotate the ECR velocity vectors for each MIF.  If there
    ! is no next record, use the difference of the last two records.
    if ( recno < nummmaf ) then
      read ( unit=in_unit, rec=recno+1 ) &
        & limb_hdr(1-w), limb_stat(1-w), limb_rad(1-w), limb_oa(1-w)
      ! swap all OA bytes:
      call swap_OA_rec ( limb_oa(1-w) )

      vel_ECI(:,1-w) = limb_oa(1-w)%sat_vel
      vNI(1-w) = norm2(vel_ECI(:,1-w))

      call sat_vel_in_ecr ( limb_oa(1-w), asciiUTC(1-w), sat_vel(:,1-w), posecr(:,1-w) )

      vN(1-w) = norm2(sat_vel(:,1-w))
      dV = acos(dot_product(sat_vel(:,w),sat_vel(:,1-w))/(vN(w)*vN(1-w))) / noMIFs
      dVN = ( vN(1-w) / vN(w) - 1.0 ) / noMIFs
      dVI = acos(dot_product(vel_ECI(:,w),vel_ECI(:,1-w))/(vNI(w)*vNI(1-w))) / noMIFs
      dVIN = ( vNI(1-w) / vNI(w) - 1.0 ) / noMIFs
write ( 42,'(a/2x,3f9.2,1p,4g15.6)' ) &
& '     VelECI   VelECR      vNI ----- dV ----- ----- dVN ---- ----- dVI ---- ---- dVIN ----', &
& vN(1-w), vNI(1-w), norm2(posecr(:,1-w))/1000, dV*rad2deg, dVN*1.0e6, dVI*rad2deg, dVIN*1.0e6
      dVV = sat_vel(:,1-w) - sat_vel(:,w)
    end if

    MAFno = MAFno + 1
write ( 42, '(1x,2(31x,a25),a)' ) AsciiUTC(w:1-w:1-2*w), &
& ' ----------------------- PosECR -----------------------'
write ( 42, '(2x,6(3f9.2:","))' ) ( limb_oa(i)%sat_vel*1000, sat_vel(:,i), i=w,1-w,1-2*w), posECR/1000

    call process_one_MAF ( MAFno, rad_sd_id, oa_sd_id, dV, dVN, dVV, &
                         & Sat_Vel(:,w:1-w:1-2*w), posECR(:,w:1-w:1-2*w), &
                         & Vel_ECI(:,w:1-w:1-2*w), dVI, dVIN, AsciiUTC, &
                         & limb_hdr(w), limb_stat(w), limb_rad(w), limb_oa(w) )

    ! Swap buffers
    w = 1 - w

  end do

  print *, 'mafs: ', MAF1, MAFn

  close (in_unit)

  call MLS_closeFile (rad_sd_id, HDFversion=HDFversion)
  call MLS_closeFile (oa_sd_id, HDFversion=HDFversion)
  call MLS_h5close (error)

9 continue

contains

  subroutine Process_One_MAF ( MAFNo, Rad_SD_id, OA_SD_id, dV, dVN, dVV, &
                             & Sat_Vel, PosECR, Vel_ECI, dVI, dVIN, AsciiUTC, &
                             & Limb_Hdr, Limb_Stat, Limb_Rad, Limb_OA )

    use OA_File_Contents, only: EMLS_OA_t
    use SDPToolkit, only: PGS_TD_UTCtoTAI
    use SwapEndian, only: SwapShort
    use Time_m, only: MS_to_HMS
    use UARS_to_EMLS_OA_m, only: UARS_to_EMLS_OA

    integer, intent(in) :: MAFNo, Rad_SD_id, OA_SD_id
    real(r8), intent(in) :: dV     ! delta V angle per MIF, ECR
    real(r8), intent(in) :: dVN    ! relative delta V length per MIF, ECR
    real(r8), intent(in) :: dVV(3) ! delta V per MAF, ECR
    real(r8), intent(in) :: Sat_Vel(3,0:1) ! Satellite velocities, ECR, km/s
    real(r8), intent(in) :: Vel_ECI(3,0:1) ! Satellite velocities, ECI, km/s
    real(r8), intent(in) :: dVI    ! delta V angle per MIF, ECI
    real(r8), intent(in) :: dVIN   ! relative delta V length per MIF, ECI
    real(r8), intent(in) :: PosECR(3,0:1)  ! Satellite positions, ECI, meters
    character(*), intent(in) :: AsciiUTC(2)
    type(limb_hdr_t), intent(inout) :: Limb_Hdr
    type(limb_stat_t), intent(inout) :: Limb_Stat
    type(limb_rad_t), intent(inout) :: Limb_Rad
    type(limb_oa_t), intent(inout) :: Limb_OA

    real :: rad(15*6,32), rad6(15,6,32), prec(15*6,32), prec6(15,6,32)
    equivalence (rad, rad6), (prec, prec6)

    integer :: DOY
    type(emls_oa_t) :: EMLS_OA
    integer :: Hrs, I, J, K, Mins, MS, Secs
    integer :: Stat ! From PGS_TD_UTCtoTAI
    integer :: Year

    character (len=*), parameter :: &
      DateFmt = "(I4, '-', I3.3, 'T', I2.2, ':', I2.2, ':', I2.2)"

    ! convert to little endian as needed:
    limb_hdr%recordno = SwapBig (limb_hdr%recordno)
    limb_hdr%mmafno = SwapBig (limb_hdr%mmafno)
    do i = 1, 2
       limb_hdr%mmaf_time(i) = SwapBig (limb_hdr%mmaf_time(i)) ! - n_days
    end do
    do i = 1, 16
       limb_stat%prd_temps(i) = SwapBig (limb_stat%prd_temps(i))
    end do
    limb_stat%maneuver_stat = SwapBig (limb_stat%maneuver_stat)
    limb_stat%mls_status = SwapBig (limb_stat%mls_status)
    do i = 1, 90
        limb_stat%wall_mmaf(i) = SwapShort (limb_stat%wall_mmaf(i))
    end do

    ! convert radiances:

    do k = 1, 32         ! MIF number
      do j = 1, 2        ! radiance/error
        do i = 1, 90     ! radiometer and channel
          limb_rad%rad_l1(i,j,k) = SwapShort (limb_rad%rad_l1(i,j,k))
        end do
      end do
    end do

    rad = limb_rad%rad_l1(:,1,:)
    prec = limb_rad%rad_l1(:,2,:)
    where ( rad < 0.0 .AND. rad /= -1000.0 ) rad = rad + 65536.0
    where ( prec < 0.0 .AND. prec /= -1000.0 ) prec = prec + 65536.0
    where ( rad /= -1000.0 ) rad = rad * 0.01   ! scale radiances

    where ( prec == -1000.0 ) prec = -1.0
    where ( prec /= -1.0 ) prec = prec * 0.01   ! scale precisions

    year = limb_hdr%mmaf_time(1) / 1000 + 1900
    doy = MOD (limb_hdr%mmaf_time(1), 1000)

    ! convert millisecs to HMS:

    call ms_to_hms ( limb_hdr%mmaf_time(2), hrs, mins, secs )
    write ( utc_time, fmt=DateFmt ) year, doy, hrs, mins, secs

    ! get and save the GIRD time:

    stat = PGS_TD_UTCtoTAI (utc_time, tai_time)
    gird_time = tai_time + TAItoGIRD

    ! write radiances etc. to file:

    call OutputL1B_rad ( rad_sd_id, mafno, limb_hdr%mmafno, rad6, prec6, &
                       & gird_time, limb_hdr )

    ! convert to EMLS OA:

    call UARS_to_EMLS_OA ( n_days, limb_oa, emls_oa, dV, dVN, dVV, Sat_Vel, &
      & posECR, Vel_ECI, dVI, dVIN, AsciiUTC )

    ! write to OA file:

    call OutputL1B_OA ( oa_sd_id, mafno, limb_hdr%mmafno, limb_oa, emls_oa )

  end subroutine Process_One_MAF

  subroutine Sat_Vel_In_ECR ( Limb_OA, AsciiUTC, Sat_vel, PosECR )
    use PGS_Interfaces, only: pgs_csc_ecitoecr, pgs_csc_ecrtoeci, pgs_csc_geotoecr

    type(limb_oa_t), intent(in) :: Limb_OA
    character(*), intent(out) :: AsciiUTC
    real(r8), intent(out) :: Sat_Vel(3)    ! ECR, m/s
    real(r8), intent(out) :: PosECR(3)     ! Meters

    real(r8) :: ECI(6,1), ECR(6,1)
    integer :: Stat

    stat = pgs_csc_geotoecr ( limb_oa%sat_long*Deg2Rad, &
                              limb_oa%sat_geod_lat*Deg2Rad, &
                              limb_oa%sat_geod_alt*1000.0d0, &
                              earthEllipseTag, posecr )

    call frame_TAI ( limb_oa, secTAI93, asciiUTC )

    ecr(1:3,1) = posECR
    ecr(4:6,1) = 0.0 ! Velocity
    ! Compute position in ECI
    stat = pgs_csc_ecrtoeci ( 1, asciiUTC, [0.0d0], ecr, eci )
    ! Now convert position and velocity to ECR
    eci(4:6,1) = limb_oa%sat_vel * 1000 ! ECI, m/s
    stat = pgs_csc_ecitoecr ( 1, asciiUTC, [0.0d0], eci, ecr )
    sat_vel = ecr(4:6,1)                ! Now in ECR, m/s

write ( 42,'(2x,2a)' ) ' ------------------------ ECI -------------------------', &
& ' ------------------------ ECR -------------------------'
write ( 42, '(2x,4(3f9.2:","))' ) eci(1:3,1)/1000,eci(4:6,1), ecr(1:3,1)/1000,ecr(4:6,1)
  end subroutine Sat_Vel_In_ECR

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end program Conv_UARS

! $Log$
! Revision 1.4  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
