! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM L2GPSpeedTest ! tests 2 ways of writing l2gp files
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use HDFEOS5, only: HE5T_NATIVE_CHAR
   use L2GPData, only: AppendL2GPData, L2GPData_T, ReadL2GPData, WriteL2GPData, &
     SetupNewL2GPRecord 
   use MLSCommon, only: R8
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
    & HDFVERSION_4, HDFVERSION_5
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStrings, only: GetStringElement, readIntsFromChars
   use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
   use PCFHdr, only: GlobalAttributes
   use Time_M, only: Time_Now, Use_Wall_Clock
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests 2 ways of writing l2gp files.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! rm -f l2gp.h4 l2gptest-4.out ; echo "1,4" | LF95.Linux/test > l2gptest-4.out
! rm -f l2gp.h5 l2gptest-5.out ; echo "1,5" | LF95.Linux/test > l2gptest-5.out

! Variables
   type (L2GPData_T)            :: l2gp, alias_l2gp
   integer                      :: hdfVersion
   integer                      :: nTimes, time, freq, level
   integer                      :: returnStatus, swfid, record_length, swid
   integer, parameter           :: nFreqs = 10
   integer, parameter           :: nLevels = 10
   character(len=*), parameter  :: l2gpFilename_h4 = '/bigdata/ahanzel/MOSS-3/L2-outputs/MLS-Aura_L2GP-H2O_V01-30-c02_1996d051.he5'
   character(len=*), parameter  :: l2gpFilename_h5 = '/bigdata/ahanzel/MOSS-3/L2-outputs/MLS-Aura_L2GP-H2O_V01-30-c02_1996d051.he5'
   character(len=*), parameter  :: newl2gpFilename = 'serialH2O.he5'
   character(len=*), parameter  :: dwl2gpFilename = 'dwH2O.he5'
   character(len=len(l2gpFilename_h4))  :: l2gpFilename
   character(len=*), parameter  :: swathName = 'H2O'
   character(len=80)            :: inputString
   character(len=8), dimension(2) :: int_chars
   integer, dimension(2)        :: the_ints
   real(r8)                     :: diff
   integer, external            :: he5_swwrlattr, he5_swopen, he5_EHwrglatt, &
     &                             HE5_SWattach, HE5_SWdetach, he5_SWsetalias, &
     &                             HE5_SWrdfld 
   real :: T0, T1, T2               ! For timing
   character(len=*), parameter  :: alphabet = &
     & 'abcdefghijklmnopqrstuvwxyz'
   character(len=40000)         :: big_file
   integer :: offset, lastProfile
   MLSMessageConfig%useToolkit = .false.
   MLSMessageConfig%LogFileUnit = -1
   ! Some extra assigments
   GlobalAttributes%InstrumentName = 'Totally bogus name'
   GlobalAttributes%ProcessLevel = 'Level 9'
   GlobalAttributes%InputVersion = 'v 2.001'
   GlobalAttributes%PGEVersion = 'v 1.1'
   GlobalAttributes%StartUTC = '01012999T00:00:00Z'
   GlobalAttributes%EndUTC = '01012999T23:59:59.9Z'
   ! Executable
   call h5open_f(returnStatus)

   hdfVersion = HDFVERSION_5
   print *, 'How many profiles/chunk (e.g., 10), hdf version (4 or 5)'
   read(*, '(a80)') inputString
   the_ints = 0
   call GetStringElement( inputString, int_chars(1), 1, .FALSE.)
   call GetStringElement( inputString, int_chars(2), 2, .FALSE.)
   call readIntsFromChars ( int_chars, the_ints, ' ,')
   nTimes = the_ints(1)
   if ( the_ints(2) /= 0 ) hdfVersion = the_ints(2)
   if ( hdfVersion == HDFVERSION_4 ) then
     l2gpFilename = l2gpFilename_h4
   else
     l2gpFilename = l2gpFilename_h5
   endif
   print *,'input string  : ', inputString
   print *,'int_chars     : ', int_chars
   print *,'ints          : ', the_ints
   print *,'hdf version   : ', hdfVersion
   print *,'file name     : ', l2gpFilename
   print *,'num profiles  : ', nTimes
   ! Now read the l2gp file
   call time_now ( t1 )
   swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &
        & record_length, DFACC_READ, FileName=l2gpFilename, &
        & hdfVersion=hdfVersion, debugOption=.false. )
   if ( returnStatus /= 0 ) then
     print *, 'Error in opening: ', trim(l2gpFilename)
     print *, 'Error number: ', returnStatus
     stop
   endif
   call ReadL2GPData ( swfid, swathName, l2gp, &
    & hdfVersion=hdfVersion )
   returnStatus = mls_io_gen_closeF('swclose', swfid, &
    & hdfVersion=hdfVersion)

   call sayTime ( 'Reading the l2gp file' )
   call time_now ( t1 )
   ! Now rewrite the l2gp file
   swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &         
    & record_length, DFACC_CREATE, FileName=newl2gpFilename, &
    & hdfVersion=hdfVersion, debugOption=.false. )
   call writeL2GPData ( l2gp, swfid, swathName, &
    & hdfVersion=hdfVersion )             
   returnStatus = mls_io_gen_closeF('swclose', swfid, &
    & hdfVersion=hdfVersion)
   call sayTime ( 'Writing the l2gp file serially' )
   ! Now do it by directWrites
   call SetupNewL2GPRecord(alias_l2gp, &
     & l2gp%nFreqs, l2gp%nLevels, nTimes)
   alias_l2gp%l2gpValue = l2gp%l2gpValue(:,:,1:nTimes)
   alias_l2gp%l2gpPrecision = l2gp%l2gpPrecision(:,:,1:nTimes)
   alias_l2gp%latitude = l2gp%latitude(1:nTimes)
   alias_l2gp%longitude = l2gp%longitude(1:nTimes)
   alias_l2gp%solarTime = l2gp%solarTime(1:nTimes)
   alias_l2gp%solarZenith = l2gp%solarZenith(1:nTimes)
   alias_l2gp%losAngle = l2gp%losAngle(1:nTimes)
   alias_l2gp%geodAngle = l2gp%geodAngle(1:nTimes)
   alias_l2gp%time = l2gp%time(1:nTimes)
   alias_l2gp%chunkNumber = l2gp%chunkNumber(1:nTimes)
   alias_l2gp%status = l2gp%status(1:nTimes)
   alias_l2gp%quality = l2gp%quality(1:nTimes)
   call time_now ( t1 )
   swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &         
    & record_length, DFACC_CREATE, FileName=dwl2gpFilename, &
    & hdfVersion=hdfVersion, debugOption=.false. )
   if ( returnStatus /= 0 ) then
     print *, 'Error in opening: ', trim(dwl2gpFilename)
     print *, 'Error number: ', returnStatus
     stop
   endif
   offset = 0
   do
     lastProfile = min(l2gp%nTimes, offset+nTimes)
     ! print *, 'Appending to swath from,to ', offset, lastProfile
     call AppendL2GPData( alias_l2gp, swfid, &
     & swathName, dwl2gpFilename, offset, hdfVersion=hdfVersion, &
     & createSwath=(offset==0) )
     offset = offset + nTimes
     if ( offset >= l2gp%nTimes) exit
   enddo
   returnStatus = mls_io_gen_closeF('swclose', swfid, &
    & hdfVersion=hdfVersion)
   call sayTime ( 'directWriting the l2gp file' )
   call h5close_f(returnStatus)
contains
  subroutine SayTime ( What )
    character(len=*), intent(in) :: What
    call time_now ( t2 )
    if ( .true. ) then
      call output ( "Total time = " )
      call output ( dble(t2), advance = 'no' )
      call blanks ( 4, advance = 'no' )
    endif
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - t1), advance = 'yes' )
  end subroutine SayTime

!==================
END PROGRAM L2GPSpeedTest
!==================

! $Log$

