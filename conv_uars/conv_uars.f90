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

  use HDF5, only:  H5GClose_f, H5GOpen_f
  use ISO_Fortran_Env, only: Output_Unit
  use MLSAuxData, ONLY: CreateGroup_MLSAuxData
  use MLSCommon, ONLY: FileNameLen
  use MLSFiles, ONLY: MLS_openFile, MLS_closeFile
  use MLSHDF5, ONLY: MakeHDF5Attribute, MLS_h5open, MLS_h5close
  use MLSKinds, only: R8
  use Output_UARS_L1B, only: OutputL1B_OA, OutputL1B_Rad
  use Rad_File_Contents, only: Lvl1_Hdr_t, Limb_Hdr_t, Limb_Stat_t, Limb_Rad_t, &
                             & Limb_OA_t, Pad_t
  use Set_Attributes_m, only: Set_Attributes
  use Swap_OA_Rec_m, only: Swap_Limb_Hdr_Rec, Swap_Limb_Stat_Rec, &
                         & Swap_Lvl1_Hdr_Rec, Swap_Rad_Rec, Swap_OA_Rec
  use UARS_Dumps, only: Show_Time
  use UARS_to_EMLS_OA_m, only: EarthEllipseTag, Frame_TAI, NoMIFs

  implicit NONE

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Is the input file name stored as an attribute in the output file?
  logical, parameter :: Input_File_Attrib = .true.

  integer, parameter :: In_Unit=101, HDFversion=5
  character (len=FileNameLen) :: Arg
  character (len=FileNameLen) :: InFile='uars.dat'
  character (len=FileNameLen) :: OA_File='OA.h5'
  character (len=FileNameLen) :: Out_Dir='/data/'
  character (len=FileNameLen) :: Rad_File='RAD.h5'
  integer :: DOY
  integer :: Error
  integer :: Grp_ID  ! of root directory of output file
  integer :: I
  integer :: InLen
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
  integer :: W ! zero or one, to select Limb_Hdr etc.
  integer :: Year
  integer :: YrDOY

  ! Inputs
  type(lvl1_hdr_t) :: Lvl1_Hdr
  type(limb_hdr_t) :: Limb_Hdr(0:1)
  type(limb_stat_t) :: Limb_Stat(0:1)
  type(limb_rad_t) :: Limb_Rad(0:1)
  type(limb_oa_t) :: Limb_OA(0:1)
  type(pad_t) :: Pad

  character (len=25) :: AsciiUTC(0:1)
  character (len=17) :: UTC_Time
  character (len=8) :: N_Str
  character (len=8) :: YrDoy_str = 'yyyydnnn'
  character (len=2048) :: CmdLine
  character (len=127) :: IOMSG
  real(r8) :: dVI, dVIN                ! delta V angle, |delta V|, ECI, per MAF
  real(r8) :: ECI(6,1), ECR(6,1)       ! Position, velocity, for conversions
  real(r8) :: GIRD_Time
  real(r8), parameter :: TAItoGIRD = 1104537627.0_r8
  real(r8) :: SecTAI93
  real(r8) :: TAI_Time
  real(r8) :: Vel_ECI(3,0:1), vNI(0:1) ! SC Velocity in ECR, its length, km/s

  ! Source code:

  ! Open UARS input file:

  MAF1 = 1
  MAFn = -1  ! Signal to get it from the header

  noargs = command_argument_count()
  if (noargs == 0) then
    print '(a)', 'No input file provided.'
    print '(2a)', 'Using defaults: Input ', trim(inFile)
    print '(2a)', 'output directory      ', trim(out_dir)
    print '(2a)', 'output OA file        ', trim(oa_file)
    print '(2a)', 'output Rad file       ', trim(rad_file)
  else
    i = 1

    do while (i <= noargs)

      call get_command_argument (i, arg)

      select case (arg)
      case default
        if ( arg(1:1) == '-' ) then
          call usage
          go to 9
        end if
        infile = TRIM(arg)                 ! input filename
        if (noargs < 3) then
          print *, 'Missing output directory argument? Too few args!'
          go to 9
        end if
      case ('-b', '--backdate')
        i = i + 1
        call get_command_argument (i, n_str)  ! how many days to backdate times in file
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
        call get_command_argument (i, n_str) ! First MAF
        read ( n_str, * ) MAF1
        i = i + 1
        call get_command_argument (i, n_str) ! Last MAF
        read ( n_str, * ) MAFn
      case ('-o', '--outdir')
        i = i + 1
        call get_command_argument (i, out_dir)  ! next argument for out dir
      case ('-h', '--help')
        call usage
        go to 9
      end select
      i = i + 1   ! Next argument
    end do
  end if

  call get_command ( cmdLine )
  print '(a)', trim(cmdLine)

  ! The length computed by
  ! inquire ( iolength=i ) limb_hdr(1), limb_stat(1), limb_rad(1), limb_oa(1)
  ! is 90 less than inlen.  What do the other 360 bytes contain?  It looks
  ! like zeroes, but is it just padding or do some records contain
  ! something useful there?

  inquire ( iolength = inlen ) &
    & limb_hdr(0), limb_stat(0), limb_rad(0), limb_oa(0), pad
  open ( unit=in_unit, file=infile, status="OLD", form="UNFORMATTED", &
       & access="DIRECT", recl=inlen, iostat=error, iomsg=iomsg )

  if (error /= 0) then
    print '(a,i0,2a/a)', 'Error status ', error, ' while opening ', trim(infile), &
      & trim(iomsg)
    go to 9
  end if

  ! read and convert header part:

  do recno = 1, 3

    read (unit=in_unit, rec=recno) lvl1_hdr

    if ( lvl1_hdr%rec_type /= 'H' ) then
      print '(a,i0,3a)', 'Something is wrong.  Record ', i, &
                       & ' is not type H.  Type is "', lvl1_hdr%rec_type, '"'
      stop
    end if
    ! convert to little endian as needed:
    call swap_lvl1_hdr_rec ( lvl1_hdr )

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
  call Swap_Limb_Hdr_Rec ( limb_hdr(0) )
  yrdoy = limb_hdr(0)%mmaf_time(1) - n_days ! This lets us backdate the files; e.g., for sids

  if ( noargs > 0 ) then   ! use passed input for setting outputs
    year = yrdoy / 1000 + 1900
    doy = mod (yrdoy, 1000)
    write (yrdoy_str, '(i4.4, "d", i3.3)') year, doy
!     print '(a,1x,i0,1x,a)', 'yrdoy str: ', yrdoy, yrdoy_str

    Rad_File = TRIM(out_dir) // '/MLS-UARS_L1BRAD_' // yrdoy_str // '.h5'
    OA_File = TRIM(out_dir) // '/MLS-UARS_L1BOA_' // yrdoy_str // '.h5'
  end if

  print '(a)',    'infile: ' // trim(infile)
  print '(a)',    'rad file: ' // trim(Rad_File)
  print '(a)',    'OA file: ' // trim(OA_File)
  print '(a,1x,i0)', 'num. days backdated:   ', n_days

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
    call h5gclose_f ( grp_id, stat )
  end if

  ! print stuff

  print '(a,1x,i0)',         'UARS day:              ', lvl1_hdr%uars_day
  print '(a,2(1x,i0),1x,a)', 'Start time:            ', lvl1_hdr%start_time, &
    & show_time ( lvl1_hdr%start_time )
  print '(a,2(1x,i0),1x,a)', 'Stop time:             ', lvl1_hdr%stop_time, &
    & show_time ( lvl1_hdr%stop_time )
  print '(a,2(1x,i0))',      'nummmaf:               ', nummmaf
  print '(a,2(1x,i0))',      'mafs:                  ', MAF1, MAFn
  ! read and convert data:

  MAFno = MAF1 - 1

  w = 0
  read ( unit=in_unit, rec=MAF1+3 ) &
    & limb_hdr(w), limb_stat(w), limb_rad(w), limb_oa(w)

  if ( limb_hdr(w)%rec_type /= 'D' ) then
    print '(a,i0,3a)', 'Something is wrong.  Record ', MAF1+3, &
                     & ' is not type D.  Type is "', limb_hdr(w)%rec_type, '"'
    stop
  end if
  ! convert to little endian as needed:
  call swap_limb_hdr_rec ( limb_hdr(w) )
  call swap_limb_stat_rec ( limb_stat(w) )
  call swap_rad_rec ( limb_rad(w) )
  call swap_OA_rec ( limb_oa(w) )

  vel_ECI(:,w) = limb_oa(w)%sat_vel
  vNI(w) = norm2(vel_ECI(:,w))

  dVI = 0.0
  dVIN = 0.0

  do recno = MAF1+3, MAFn+3

    ! Get the record for the next MAF, if there is one.
    ! Calculate the change in velocity vector angle.
    ! This is used to rotate the ECR velocity vectors for each MIF.  If there
    ! is no next record, use the difference of the last two records.
    if ( recno < nummmaf+3 ) then
      read ( unit=in_unit, rec=recno+1 ) &
        & limb_hdr(1-w), limb_stat(1-w), limb_rad(1-w), limb_oa(1-w)
      if ( limb_hdr(1-w)%rec_type /= 'D' ) then
        print '(a,i0,3a)', 'Something is wrong.  Record ', recno, &
                         & ' is not type D.  Type is "', limb_hdr(1-w)%rec_type, '"'
        stop
      end if
      ! convert to little endian as needed:
      call swap_limb_hdr_rec ( limb_hdr(1-w) )
      call swap_limb_stat_rec ( limb_stat(1-w) )
      call swap_rad_rec ( limb_rad(1-w) )
      call swap_OA_rec ( limb_oa(1-w) )

      vel_ECI(:,1-w) = limb_oa(1-w)%sat_vel
      vNI(1-w) = norm2(vel_ECI(:,1-w))

      dVI = acos(dot_product(vel_ECI(:,w),vel_ECI(:,1-w))/(vNI(w)*vNI(1-w))) / noMIFs
      dVIN = ( vNI(1-w) / vNI(w) - 1.0 ) / noMIFs
    end if

    MAFno = MAFno + 1

    call process_one_MAF ( MAFno, rad_sd_id, oa_sd_id, &
                         & dVI, dVIN, &
                         & limb_hdr(w), limb_stat(w), limb_rad(w), limb_oa(w) )

    ! Swap buffers
    w = 1 - w

  end do

  close (in_unit)

  call MLS_closeFile (rad_sd_id, HDFversion=HDFversion)
  call MLS_closeFile (oa_sd_id, HDFversion=HDFversion)
  call MLS_h5close (error)

9 continue

contains

  subroutine Process_One_MAF ( MAFNo, Rad_SD_id, OA_SD_id, &
                             & dVI, dVIN, &
                             & Limb_Hdr, Limb_Stat, Limb_Rad, Limb_OA )

    use OA_File_Contents, only: EMLS_OA_t
    use SDPToolkit, only: PGS_SMF_GetMsg, PGS_S_Success, PGS_TD_UTCtoTAI
    use Time_m, only: MS_to_HMS
    use UARS_to_EMLS_OA_m, only: UARS_to_EMLS_OA

    integer, intent(in) :: MAFNo, Rad_SD_id, OA_SD_id
    real(r8), intent(in) :: dVI    ! delta V angle per MIF, ECI
    real(r8), intent(in) :: dVIN   ! relative delta V length per MIF, ECI
    type(limb_hdr_t), intent(inout) :: Limb_Hdr
    type(limb_stat_t), intent(inout) :: Limb_Stat
    type(limb_rad_t), intent(inout) :: Limb_Rad
    type(limb_oa_t), intent(inout) :: Limb_OA
    character(255) :: Err, Msg

    real :: rad(15*6,32), rad6(15,6,32), prec(15*6,32), prec6(15,6,32)
    equivalence (rad, rad6), (prec, prec6)

    integer :: DOY
    type(emls_oa_t) :: EMLS_OA
    integer :: Hrs, I, J, K, Mins, MS, Secs
    integer :: Stat ! From PGS_TD_UTCtoTAI
    integer :: Year

    character (len=*), parameter :: &
      DateFmt = "(I4, '-', I3.3, 'T', I2.2, ':', I2.2, ':', I2.2)"

    ! convert radiances:

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
    if ( stat /= PGS_S_Success ) then
      print '(a,i0)', 'Status from PGS_TD_UTCtoTAI = ', stat
      print '(a)', 'Is the environment variable PGSHOME set?'
      call PGS_SMF_GetMsg ( stat, err, msg )
      print '(a)', trim(err), trim(msg)
      stop
    end if
    gird_time = tai_time + TAItoGIRD

    ! write radiances etc. to file:

    call OutputL1B_rad ( rad_sd_id, mafno, limb_hdr%mmafno, rad6, prec6, &
                       & gird_time, limb_hdr )

    ! convert to EMLS OA:

    call UARS_to_EMLS_OA ( n_days, limb_oa, emls_oa, dVI, dVIN )

    ! write to OA file:

    call OutputL1B_OA ( oa_sd_id, mafno, limb_hdr%mmafno, limb_oa, emls_oa )

  end subroutine Process_One_MAF

  subroutine Usage
    character(255) :: Line
    call get_command_argument ( 0, line )
    print '(3a)', 'Usage: ', trim(line), ' [options] [input [output directory]]'
    print '( a)', ' Options: -b n_days           => number of days to backdate'
    print '( a)', '          --backdate n_days   => number of days to backdate'
    print '( a)', '          -M first last       => MAFs to run'
    print '( a)', '          --MAF first last    => MAFs to run'
    print '( a)', '          -o output_dir       => output directory'
    print '( a)', '          --outdir output_dir => output directory'
    print '( a)', '          -h                  => print help'
    print '( a)', '          --help              => print help'
    print '( a)', ' Defaults: options: None'
    print '(2a)', '  input file:       ', trim(inFile)
    print '(2a)', '  output directory: ', trim(out_dir)
    print '(2a)', '  OA file:          ', trim(oa_file)
    print '(2a)', '  RAD file:         ', trim(rad_file)
  end subroutine Usage

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
! Revision 1.9  2016/12/21 20:37:41  vsnyder
! Fiddle some command-line stuff
!
! Revision 1.8  2015/04/21 02:27:08  vsnyder
! Create attributes for some components' units
!
! Revision 1.7  2015/04/21 01:18:54  vsnyder
! Use swap routines in Swap_OA_Rec_m.  Get record size using INQUIRE with
! IOLENGTH=.  Spiff some printing.  Print usage with -h or --help option.
! Check record types.  Calculate highest record number properly.
!
! Revision 1.6  2015/01/24 02:09:16  vsnyder
! MIF resolve most quantities
!
! Revision 1.5  2015/01/22 02:18:56  vsnyder
! PGS_Interfaces.f90
!
! Revision 1.4  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
