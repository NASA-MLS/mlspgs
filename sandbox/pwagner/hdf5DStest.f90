! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM hdf5DStest ! tests MLSHDF5 routines
!=================================

   use dump_0, only: dump
   use Hdf, only: DFACC_CREATE, DFACC_RDWR, DFACC_READ
   use MLSCommon, only: r4
   USE MLSFILES, ONLY: MLS_SFSTART, MLS_SFEND, &
    &  HDFVERSION_4, HDFVERSION_5
   use MLSHDF5, only: SaveAsHDF5DS, LoadFromHDF5DS, MakeHDF5Attribute, &
    &  GetHDF5Attribute
   use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMessageConfig
   use MLSStrings, only: GetStringElement, readIntsFromChars
   use HDF5, only: h5gclose_f, h5gopen_f
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the l2auxData subroutines.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! rm -f l2aux.h4 l2auxtest.out ; echo "1,4" | LF95.Linux.atlas/test > l2auxtest.out

! Variables
   integer                      :: hdfVersion
   integer                      :: nTimes, time, burst, level, burst_length
   integer                      :: returnStatus, fileID, grp_id
   integer                      :: end_time, start_time
   integer, parameter           :: nBursts = 10
   integer, parameter           :: nLevels = 10
   character(len=*), parameter  :: l2auxFilename_h4 = 'l2aux.h4'
   character(len=*), parameter  :: l2auxFilename_h5 = 'l2aux.h5'
   character(len=len(l2auxFilename_h4))  :: l2auxFilename
   character(len=*), parameter  :: swathName = 'test swath'
   character(len=*), parameter  :: GlAttrN = 'global attribute'
   character(len=*), parameter  :: GlAttrV = 'GA value'
   character(len=*), parameter  :: DSname = 'test data'
   character(len=*), parameter  :: DSUnits = 'test units'
   character(len=*), parameter  :: DimName = 'test dimension'
   character(len=*), parameter  :: DimUnits = 'dimension units'
   real(r4), dimension(nLevels) :: DimValues
   character(len=80)            :: inputString
   character(len=8), dimension(2) :: int_chars
   integer, dimension(2)        :: the_ints
   real(r4)                     :: diff
   real(r4), dimension(:,:), pointer   :: l2auxValue => NULL()
   logical, parameter           :: WRITEATTRIBUTESDURING = .TRUE.
   ! Executable
   call h5open_f(returnStatus)
   MLSMessageConfig%useToolkit = .false.
   MLSMessageConfig%LogFileUnit = -1
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
     l2auxFilename = l2auxFilename_h4
   else
     l2auxFilename = l2auxFilename_h5
   endif
   print *,'input string  : ', inputString
   print *,'int_chars     : ', int_chars
   print *,'ints          : ', the_ints
   print *,'hdf version   : ', hdfVersion
   print *,'file name     : ', l2auxFilename
   print *,'num profiles  : ', nTimes
   do level=1, nLevels
     DimValues(level) = (level-.5) / nLevels
   enddo
!  Loop of bursts
   end_time = 0
   burst_length = 1 + (nTimes-1)/nBursts
   do burst=1, nBursts
     start_time = end_time + 1
     end_time = min(end_time + burst_length, nTimes)
     if ( end_time < start_time) exit
     print *,'start_time    : ', start_time
     print *,'end_time      : ', end_time
     allocate(l2auxValue(nLevels, 1:end_time-start_time+1), stat=returnStatus)
     if ( returnStatus /= 0 ) then
       print *, 'Error allocating l2auxValue'
       stop
     endif
     do time=start_time, end_time
       do level=1, nLevels
           l2auxValue(level, time-start_time+1) = time + 100*level
       enddo
     enddo
     if ( burst == 1 ) then
       fileID = mls_sfstart ( trim(l2auxFilename), DFACC_CREATE, &
         &                                          hdfVersion=hdfVersion )
     else
       fileID = mls_sfstart ( trim(l2auxFilename), DFACC_RDWR, &
         &                                          hdfVersion=hdfVersion )
     endif
     call SaveAsHDF5DS ( fileID, DSname, l2auxValue, &
      & (/0, start_time-1/), (/nLevels, end_time-start_time+1/), &
      & may_add_to=.true., adding_to=(burst /= 1) )             

     if ( WRITEATTRIBUTESDURING ) then
       call MakeHDF5Attribute ( fileID, DSName, 'Units', DSUnits, &
        & skip_if_already_there=.true. )
       call MakeHDF5Attribute ( fileID, DSName, 'Dimension', DimName, &
        & skip_if_already_there=.true. )
       call MakeHDF5Attribute ( fileID, DSName, 'DimesnionValues', DimValues, &
        & skip_if_already_there=.true. )
       call MakeHDF5Attribute ( fileID, DSName, 'DimesnionUnits', DimUnits, &
        & skip_if_already_there=.true. )
     endif

     returnStatus = mls_sfend( fileID, hdfVersion=hdfVersion)
     if ( returnStatus /= 0 ) then
       print *, 'Error in ending hdf access to file'
       stop
     endif
     deallocate(l2auxValue, stat=returnStatus)
     if ( returnStatus /= 0 ) then
       print *, 'Error deallocating l2auxValue'
       stop
     endif
   enddo
   ! Now read what we just wrote and see if it's the same
   fileID = mls_sfstart ( trim(l2auxFilename), DFACC_READ, &
         &                                          hdfVersion=hdfVersion )
   allocate(l2auxValue(nLevels, nTimes), stat=returnStatus)
   if ( returnStatus /= 0 ) then
     print *, 'Error allocating l2auxValue'
     stop
   endif
   call LoadFromHDF5DS (fileID, DSname, l2auxValue)
   returnStatus = mls_sfend( fileID, hdfVersion=hdfVersion)

   ! Dump 2 l2auxs and accumulate differences
   call dump(l2auxValue, 'test data')
   deallocate(l2auxValue, stat=returnStatus)
   if ( returnStatus /= 0 ) then
     print *, 'Error deallocating l2auxValue'
     stop
   endif
   
   ! Now write attributes to data
   fileID = mls_sfstart ( trim(l2auxFilename), DFACC_RDWR, &
     &                                          hdfVersion=hdfVersion )
   if ( .not. WRITEATTRIBUTESDURING ) then
     call MakeHDF5Attribute ( fileID, DSName, 'Units', DSUnits )
     call MakeHDF5Attribute ( fileID, DSName, 'Dimension', DimName )
     call MakeHDF5Attribute ( fileID, DSName, 'DimesnionValues', DimValues )
     call MakeHDF5Attribute ( fileID, DSName, 'DimesnionUnits', DimUnits )
   endif
   ! Now try to write global attributes
   call h5gopen_f(fileID, '/', grp_id, returnStatus)
   call MakeHDF5Attribute ( grp_id, GlAttrN, GlAttrV)
   call h5gclose_f(grp_id, returnStatus)
   returnStatus = mls_sfend( fileID, hdfVersion=hdfVersion)

   call h5close_f(returnStatus)
!==================
END PROGRAM hdf5DStest
!==================

! $Log$

