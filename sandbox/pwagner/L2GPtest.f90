! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM L2GPtest ! tests L2GPData routines
!=================================

   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use L2GPData, only: Dump, L2GPData_T, ReadL2GPData, WriteL2GPData, &
     SetupNewL2GPRecord 
   use MLSCommon, only: R8
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
    & HDFVERSION_4, HDFVERSION_5
   use MLSStrings, only: GetStringElement, readIntsFromChars
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the L2GPData subroutines.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! rm -f l2gp.h4 l2gptest.out ; echo "1,4" | LF95.Linux.atlas/test > l2gptest.out

! Variables
   type (L2GPData_T)            :: l2gp, read_l2gp
   integer                      :: hdfVersion
   integer                      :: nTimes, time, freq, level
   integer                      :: returnStatus, swfid, record_length
   integer, parameter           :: nFreqs = 10
   integer, parameter           :: nLevels = 10
   character(len=*), parameter  :: l2gpFilename_h4 = 'l2gp.h4'
   character(len=*), parameter  :: l2gpFilename_h5 = 'l2gp.h5'
   character(len=len(l2gpFilename_h4))  :: l2gpFilename
   character(len=*), parameter  :: swathName = 'test swath'
   character(len=80)            :: inputString
   character(len=8), dimension(2) :: int_chars
   integer, dimension(2)        :: the_ints
   real(r8)                     :: diff
   ! Executable
   call h5open_f(returnStatus)

   hdfVersion = HDFVERSION_5
   print *, 'How many profiles (e.g., 100), hdf version (4 or 5)'
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
   call SetupNewL2GPRecord(l2gp, nFreqs, nLevels, nTimes)
   l2gp%name = swathName
   do time=1, nTimes
     do level=1, nLevels
       do freq=1, nFreqs
         l2gp%l2gpValue(freq, level, time) = freq + 100*level
       enddo
     enddo
     l2gp%latitude(time) = time
     l2gp%longitude(time) = time
     l2gp%solarTime(time) = time
     l2gp%solarZenith(time) = time
     l2gp%losAngle(time) = time
     l2gp%geodAngle(time) = time
     l2gp%time(time) = time
     l2gp%chunkNumber(time) = time
   enddo
   l2gp%l2gpPrecision = l2gp%l2gpValue/10000
   do level=1, nLevels
     l2gp%pressures(level) = level
   enddo
   do freq=1, nFreqs
     l2gp%frequency(freq) = freq
   enddo
   print *, 'Frequencies: (before writeL2GPData)'
   print *, l2gp%frequency
   swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &         
    & record_length, DFACC_CREATE, FileName=l2gpFilename, &
    & hdfVersion=hdfVersion, debugOption=.false. )

   call writeL2GPData ( l2gp, swfid, swathName, &
    & hdfVersion=hdfVersion )             

   print *, 'Frequencies: (after writeL2GPData)'
   print *, l2gp%frequency
   returnStatus = mls_io_gen_closeF('swclose', swfid, &
    & hdfVersion=hdfVersion)
   ! Now read what we just wrote and see if it's the same
   swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &
        & record_length, DFACC_READ, FileName=l2gpFilename, &
        & hdfVersion=hdfVersion, debugOption=.false. )
   call ReadL2GPData ( swfid, swathName, read_l2gp, &
    & hdfVersion=hdfVersion )
   read_l2gp%name = swathName // ' (read)'
   returnStatus = mls_io_gen_closeF('swclose', swfid, &
    & hdfVersion=hdfVersion)

   ! Dump 2 l2gps and accumulate differences
   call dump(l2gp)
   call dump(read_l2gp)
   diff = 0.
   do time=1, nTimes
     do level=1, nLevels
       do freq=1, nFreqs
         diff = diff + &
          & abs( &
          & l2gp%l2gpValue(freq, level, time) &
          & - &
          & read_l2gp%l2gpValue(freq, level, time) &
          & )
       enddo
     enddo
   enddo
   print *, 'diff: ', diff
   call h5close_f(returnStatus)
!==================
END PROGRAM L2GPtest
!==================

! $Log$
! Revision 1.1  2003/01/15 19:25:56  pwagner
! First commit
!

