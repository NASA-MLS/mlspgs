program sids_l1boa ! Generate simulated L1BOA data for MLS

  use Hdf, only: DFACC_WRITE, SFEND, SFSTART
  use L1boa, only: L1BOA_FILL, L1BOA_NOFILL
  use MLSCommon, only: R8, TAI93_RANGE_T
  use MLSMessageModule, only: MLSMSG_ALLOCATE, MLSMSG_ERROR, MLSMSG_FILEOPEN, &
    & MLSMESSAGE
  use Orbit, only: ALTG, ALTT, ASCTAI, DSCTAI, ORBINCLINE, ORBITNUMBER, &
    & ORBIT_MET, NUMORB
  use OutputL1B, only: L1BOAINDEX_T
  use Read, only: READ_SCAN, READ_UIF
  use Sd, only: SD_CREATE
  use SDPToolkit, only: PGS_PC_GETREFERENCE
  use Time, only: TIME_DAY, TIME_MAF, TIME_POST, TIME_PRE

  implicit none

  !------------------- RCS Ident Info -----------------------
  character(len=130) :: Id = &                                                    
    "$Id$"
  character(len=*), parameter :: ModuleName= "$RCSfile$"
  !----------------------------------------------------------

  ! Brief description of program
  ! This program creates a  SIDS L1BOA file.

  ! Parameters
  integer, parameter :: SD_OUTFILE = 333

  ! Variables
  type ( TAI93_Range_T ) :: times
  type( L1BOAindex_T) :: index

  character(len=27) :: endTime, startTime
  character(len=135) :: physicalFilename
  character(len=480) :: msr
  character(len=27), allocatable :: mafTime(:), postTimes(:), preTimes(:)
  character(len=27), allocatable :: scanTimes(:,:)

  integer :: alloc_err, dealloc_err, i, lenMAF, lenMIF, mifG, mifRate, mifT
  integer :: minSC, nPhaseG, nPhaseT, offsetMAF, opd, orbPerDay
  integer :: postMAF, preMAF, returnStatus, scansPerOrb, sd_id, status
  integer :: SD_VERSION
  integer, allocatable :: numMIFs(:), nV(:)
  integer, allocatable :: numValues(:,:)

  real, allocatable :: scRate(:), scRateT(:)

  real(r8) :: homeAlt, homeLat, mfps, spm
  real(r8), allocatable :: mafTAI(:), offsets(:), postTAI(:), preTAI(:)
  real(r8), allocatable :: scanTAI(:,:)

  ! Read SIDS input
  call Read_uif(altG, altT, times, endTime, homeAlt, homeLat, mifG, mifT, &
    mifRate, nPhaseG, nPhaseT, offsetMAF, scansPerOrb, startTime)

  allocate(scRate(mifG), scRateT(mifT), STAT=alloc_err)
  if ( alloc_err /= 0 ) then
    msr = MLSMSG_Allocate // '  scan rate variables.'
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  call Read_scan(mifG, mifT, nPhaseG, nPhaseT, scRate, scRateT)

  ! Get orbit metadata for entire day
  call Orbit_met(startTime, times, ascTAI, dscTAI, numOrb, orbitNumber)

  opd = numOrb+1

  ! Calculate MAF start times
  allocate( numMIFs(scansPerOrb), numValues(scansPerOrb,opd), &
    postTAI(scansPerOrb), postTimes(scansPerOrb), preTAI(scansPerOrb), &
    preTimes(scansPerOrb), scanTAI(scansPerOrb,opd), &
    scanTimes(scansPerOrb,opd), STAT=alloc_err )
  if ( alloc_err /= 0 ) then
    msr = MLSMSG_Allocate // '  MAF quantities.'
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  ! Need to zero numValues. The program does not set all the elements of 
  ! numValues(1:orbPerDay-1) before it calculates minSC. This may lead to 
  ! minSC being (e.g.) -28484737282 yucky consequences. This is doubtless
  ! compiler dependent -- it is needed for NAG f95 / linux. (HCP)
  numValues=0

  mfps = dble(mifRate)
  spm = 1/mfps

  call Time_MAF(startTime, times, homeAlt, homeLat, mifRate, scansPerOrb, &
    spm, scanTimes, scanTAI, numValues, orbPerDay)

  call Time_pre(numValues(:,1), scanTimes(1,1), times%startTime, scansPerOrb, &
    spm, preTimes, preTAI, numMIFs, preMAF)

  call Time_post(numValues(:,orbPerDay-1), scanTimes(1,orbPerDay), &
    times%endTime, scansPerOrb, spm, postTimes, postTAI, postMAF)

  ! Check the lengths of the given scans against the minimum # of MIFs per MAF
  minSC = minval( numValues(:,1:orbPerDay-1) )

  if ( mifG > minSC ) then
    print *, 'Desired GHz scan program is too long:'
    print *, 'Minimum number of MIFs per MAF = ', minSC
    print *, 'Number of MIFs in GHz scan = ', mifG
    stop
  endif

  if ( mifT > minSC ) then
    print *, 'Desired THz scan program is too long:'
    print *, 'Minimum number of MIFs per MAF = ', minSC
    print *, 'Number of MIFs in THz scan = ', mifT
    stop
  endif

  ! Calculate one-time values, outside of MAF/orbit loop(s)
  lenMIF = maxval(numValues)

  allocate( offsets(lenMIF), STAT=alloc_err )
  if ( alloc_err /= 0 ) then
    msr = MLSMSG_Allocate // '  offsets.'
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  do i = 1, lenMIF
    offsets(i) = (i-1)/mfps
  enddo

  orbIncline = 98.145

  ! Get filename from PCF
  SD_VERSION = 1

  returnStatus = Pgs_pc_getReference(SD_OUTFILE, SD_VERSION, physicalFilename)

  ! Open HDF-SD output file and create empty arrays; calculate # of MAFs in the
  ! day.
  call Sd_create(mifG, lenMIF, mifT, offsetMAF, orbPerDay, physicalFilename, &
    postMAF, preMAF, scansPerOrb, lenMAF)

  print *, 'Num of MAFs in file = ', lenMAF

  ! After creation, re-open the HDF file and initialize the SD interface for
  ! writing
  sd_id = sfstart(physicalFilename, DFACC_WRITE)
  if (sd_id == -1) then
    msr = MLSMSG_Fileopen // physicalFilename
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  ! Re-arrange time, number of values into single array for entire file
  allocate( mafTAI(lenMAF), mafTime(lenMAF), nV(lenMAF), STAT=alloc_err )
  if ( alloc_err /= 0 ) then
    msr = MLSMSG_Allocate // '  time variables.'
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  call Time_day(lenMAF, numMIFs, numValues, orbPerDay, postMAF, postTAI, &
    postTimes, preMAF, preTAI, preTimes, scansPerOrb, scanTAI, &
    scanTimes, mafTAI, mafTime, nV)

  ! For first MAF, calculate & write L1BOA data with fill values/settings
  print *, 'Processing MAF 1 ... '

  index%MAFStartTimeUTC = mafTime(1)
  index%MAFStartTimeTAI = mafTAI(1)
  index%noMIFs = nV(1)
  index%counterMAF = 0

  call L1boa_fill(altG, altT, ascTAI, dscTAI, sd_id, mifG, mifT, index, 1, &
    numOrb, offsets, orbIncline, orbitNumber, scRate, scRateT)

  ! Close the file to set fill values
  status = sfend(sd_id)
  if (status == -1) then
    msr = 'Failed to close file ' // physicalFilename
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  ! Re-open the HDF file and initialize the SD interface for writing
  sd_id = sfstart(physicalFilename, DFACC_WRITE)
  if (sd_id == -1) then
    msr = MLSMSG_Fileopen // physicalFilename
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  ! For remaining MAFs, calculate & write L1BOA data
  do i = 2, lenMAF

    print *, 'Processing MAF ', i, ' ... '

    index%MAFStartTimeUTC = mafTime(i)
    index%MAFStartTimeTAI = mafTAI(i)
    index%noMIFs = nV(i)
    index%counterMAF = i-1

    call L1boa_nofill(altG, altT, ascTAI, dscTAI, sd_id, mifG, mifT, index, &
      i, numOrb, offsets, orbIncline, orbitNumber, scRate, scRateT)

  enddo

  ! Terminate access to the SD interface and close the file
  status = sfend(sd_id)
  if (status == -1) then
    msr = 'Failed to close file ' // physicalFilename
    call MLSMessage(MLSMSG_Error, ModuleName, msr)
  endif

  deallocate(offsets, numMIFs, numValues, postTAI, postTimes, preTAI, &
    preTimes, scanTAI, scanTimes, scRate, scRateT, mafTAI, mafTime, &
    nV, STAT=dealloc_err )
  if ( dealloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName,&
    & 'Failed deallocation.')

  ! Detailed description of program
  ! The program creates a SIDS L1BOA SDS-HDF file.  It uses routines such that
  ! the fill value/modes are set only on writing the first record.

end program sids_l1boa

! $Log$
! Revision 1.6  2001/10/11 23:27:23  livesey
! Tried to change the chmod stuff
!
! Revision 1.5  2001/10/11 23:25:19  livesey
! Tidied up a bit
!
! Revision 1.4  2001/10/11 22:07:59  livesey
! Just tidied up a bit
!
! Revision 1.3  2001/10/10 23:20:12  pwagner
! Now compiles with latest l1/Orbit
!
! Revision 1.2  2000/12/04 18:36:36  pumphrey
! fixed a tiny bug
!
! Revision 1.1  2000/11/30 16:34:27  nakamura
! The SIDS L1BOA program.
