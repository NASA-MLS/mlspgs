
! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module read

  use MLSCommon, only: R8, TAI93_RANGE_T
  use SDPToolkit, only: PGS_S_SUCCESS, PGSD_IO_GEN_RSEQFRM
  implicit none
  private

  public :: READ_UIF, READ_SCAN

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &                                                    
    "$Id$"
  !----------------------------------------------------------
  ! This module contains subroutines related to reading the SIDS L1BOA
  ! User Input File.

  ! Parameters
  character (LEN=*), parameter :: labelFmt = "(A77)"
  integer, parameter :: SIDS_UIF = 322

contains

  !------------------------------------------------Read_UIF -----
  subroutine Read_uif(altG, altT, times, endTime, homeAlt, homeLat, mifG, &
    mifT, mifRate, nPhaseG, nPhaseT, offsetMAF, &
    scansPerOrb, startTime)
    ! This subroutine reads the first (non-scan) portion of the UIF and skims the
    ! scan portion for information needed to allocate variables.
    
    ! Arguments
    type( TAI93_Range_T ), intent(OUT) :: times
    character (LEN=27), intent(OUT) :: endTime, startTime
    integer, intent(OUT) :: mifG, mifT, nPhaseG, nPhaseT
    integer, intent(OUT) :: mifRate, offsetMAF, scansPerOrb
    real(r8), intent(OUT) :: altG, altT, homeAlt, homeLat

    ! Parameters
    character (LEN=*), parameter :: dtFmt = "(A10, 1X, A15)"

    ! Functions
    integer :: Pgs_io_gen_openF, Pgs_io_gen_closeF, Pgs_td_utcToTAI

    ! Variables
    character (LEN=10) :: startDate, endDate
    character (LEN=15) :: startUTC, endUTC
    character (LEN=32) :: mnemonic
    character (LEN=77) :: label
    character (LEN=480) :: msg
    integer :: dur, i, ios, j, processUIF, returnStatus, version

    ! Executable code

    version = 1

    ! Open the UIF as a generic file for reading
    returnStatus = Pgs_io_gen_openF (SIDS_UIF, PGSd_IO_Gen_RSeqFrm, 0, &
      processUIF, version)
    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      print *, 'Error opening UIF:  ', mnemonic
      print *, msg
    endif

    ! Read date/time input
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    if (ios /= 0) print *, 'Error reading UIF file.'
    read(UNIT=processUIF, IOSTAT=ios, FMT=dtFmt) startDate, startUTC

    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=dtFmt) endDate, endUTC

    ! Convert to asciiUTC, TAI formats
    startTime = startDate // 'T' // startUTC // 'Z'
    endTime = endDate // 'T' // endUTC // 'Z'

    returnStatus = Pgs_td_utcToTAI(startTime, times%startTime)
    returnStatus = Pgs_td_utcToTAI(endTime, times%endTime)

    ! Read rest of non-scan input
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) mifRate
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) scansPerOrb
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) offsetMAF
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) homeLat
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) homeAlt

    ! Read GHz scan program input
    do i = 1, 4
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      if ( label(17:24) == 'Geod Alt' ) exit
    enddo

    read(UNIT=processUIF, IOSTAT=ios, FMT=*) altG
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) nPhaseG

    ! Calculate GHz scan length
    mifG = 0
    do i = 1, nPhaseG
      do j = 1, 4
        read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
        if ( label(9:11) == 'Dur' ) exit
      enddo
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) dur
      mifG = mifG + dur
    enddo

    ! Read THz scan program input
    do i = 1, 4
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      if ( label(17:24) == 'Geod Alt' ) exit
    enddo

    read(UNIT=processUIF, IOSTAT=ios, FMT=*) altT

    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
    read(UNIT=processUIF, IOSTAT=ios, FMT=*) nPhaseT

    ! Calculate THz scan length
    mifT = 0
    do i = 1, nPhaseT
      do j = 1, 4
        read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
        if ( label(9:11) == 'Dur' ) exit
      enddo
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) dur
      mifT = mifT + dur
    enddo

    ! Close UIF
    returnStatus = Pgs_io_gen_closeF (processUIF)
    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      print *, 'Error closing UIF:  ', mnemonic
      print *, msg
    endif

    ! Convert km (input) to m (used by Toolkit) for altitudes

    homeAlt = homeAlt * 1000
    altG = altG * 1000
    altT = altT * 1000

  end subroutine Read_uif

  !------------------------------------------- Read_Scan ---------------------
  subroutine Read_scan(mifG, mifT, nPhaseG, nPhaseT, scanRate, scanRateT)
    ! This subroutine reads the second (scan) portion of the UIF and calculates the
    ! GHz and THz scan rates.

    ! Arguments
    integer, intent(IN) :: mifG, mifT, nPhaseG, nPhaseT

    real, intent(OUT) :: scanRate(mifG), scanRateT(mifT)

    ! Functions
    integer :: Pgs_io_gen_openF, Pgs_io_gen_closeF

    ! Variables
    character (LEN=32) :: mnemonic
    character (LEN=77) :: label
    character (LEN=480) :: msg
    integer :: i, ios, processUIF, returnStatus, sum, version
    integer :: dur(nPhaseG)
    integer :: durT(nPhaseT)
    real(r8) :: rate(nPhaseG)
    real(r8) :: rateT(nPhaseT)

    ! Executable code

    version = 1

    ! Open the UIF as a generic file for reading
    returnStatus = Pgs_io_gen_openF (SIDS_UIF, PGSd_IO_Gen_RSeqFrm, 0, &
      processUIF, version)
    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      print *, 'Error opening UIF:  ', mnemonic
      print *, msg
    endif

    ! Find scan program phase section
    do i = 1, 20
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      if ( label(16:22) == 'Phases' ) exit
    enddo

    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label

    ! Read GHz scan information
    sum = 0
    do i = 1, nPhaseG
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) rate(i)
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) dur(i)

      ! Calculate GHz rate information
      scanRate(sum+1:sum+dur(i)) = rate(i)
      sum = sum + dur(i)
    enddo

    ! Find scan program phase section for THz
    do i = 1, 5
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      if ( label(16:22) == 'Phases' ) exit
    enddo
    read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label

    ! Read THz scan information
    sum = 0
    do i = 1, nPhaseT
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) rateT(i)
      read(UNIT=processUIF, IOSTAT=ios, FMT=labelFmt) label
      read(UNIT=processUIF, IOSTAT=ios, FMT=*) durT(i)
      ! Calculate THz rate information
      scanRateT(sum+1:sum+durT(i)) = rateT(i)
      sum = sum + durT(i)
    enddo

    ! Close UIF
    returnStatus = Pgs_io_gen_closeF (processUIF)
    if (returnStatus /= PGS_S_SUCCESS) then
      call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
      print *, 'Error closing UIF:  ', mnemonic
      print *, msg
    endif

  end subroutine Read_scan

end module read

! $Log$
! Revision 1.1  2000/11/30 16:28:18  nakamura
! Module for reading the SIDS L1BOA User Input File.
