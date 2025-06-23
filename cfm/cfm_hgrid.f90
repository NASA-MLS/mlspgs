! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_HGrid_m
   use HGridsDatabase, only: HGrid_T, CreateEmptyHGrid, TrimHGrid, &
                             FindClosestMatch
   use MLSKinds, only: R8
   use L1BData, only: L1BData_T, ReadL1BData, DeallocateL1BData, AssembleL1BQtyName
   use MLSCommon, only: NameLen, MLSFile_T
   use Constants, only: Deg2Rad, Rad2Deg
   use MLSNumerics, only: Hunt, InterpolateValues
   use EmpiricalGeometry, only: EmpiricalLongitude, ChooseOptimumLon0
   use Allocate_Deallocate, only : allocate_test, deallocate_test
   use MLSFiles, only: GetMLSFileByType
   use Chunks_m, only: MLSChunk_T
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error

   implicit none

   public :: CreateRegularHGrid

   interface CreateRegularHGrid
       module procedure CreateRegularHGrid_maf, CreateRegularHGrid_chunk
   end interface

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   private

   contains

   ! Use the input given, and the L1BOA file to create a regular hGrid.
   type(HGrid_T) function CreateRegularHGrid_maf (instrumentModuleName, origin, &
                  spacing, insetOverLaps, filedatabase, &
                  firstL1Maf, lastL1Maf) result(hgrid)

       ! The (x,y) origin of coordinate
      real(r8), intent(in) :: origin
      real(r8), intent(in) :: spacing
      ! the filedatabase must contains L1BOA file (see CFM_MLSSetup)
      type(MLSFile_T), dimension(:), pointer :: filedatabase
      ! Either "GHz" or "THz" but the string is case-insensitive.
      character(len=*), intent(in) :: instrumentModuleName
      logical, intent(in) :: insetOverlaps
      ! Start of the maf range in the given L1BOA in the filedatabase
      integer, intent(in) :: firstL1Maf
      ! End of the maf range in the given L1BOA in the filedatabase
      integer, intent(in) :: lastL1Maf

      ! Contains the time range of the data to be gotten out of L1BOA
      ! (see CFM_MLSSetup)
      type(MLSChunk_T) :: chunk

      chunk%firstMafIndex = firstL1Maf
      chunk%lastMafIndex = lastL1Maf

      hgrid = CreateRegularHGrid_chunk(instrumentModuleName, origin, &
              spacing, insetOverlaps, filedatabase, chunk)
   end function

   ! Use the input given, and the L1BOA file to create a regular hGrid.
   type(HGrid_T) function CreateRegularHGrid_chunk (instrumentModuleName, origin, &
                  spacing, insetOverLaps, filedatabase, &
                  chunk) result(hgrid)
      ! The (x,y) origin of coordinate
      real(r8), intent(in) :: origin
      real(r8), intent(in) :: spacing
      ! the filedatabase must contains L1BOA file (see CFM_MLSSetup)
      type(MLSFile_T), dimension(:), pointer :: filedatabase
      ! Either "GHz" or "THz" but the string is case-insensitive.
      character(len=*), intent(in) :: instrumentModuleName
      logical, intent(in) :: insetOverlaps
      ! Contains the time range of the data to be gotten out of L1BOA
      ! (see CFM_MLSSetup)
      type(MLSChunk_T) :: chunk

      type(MLSFile_T), pointer :: l1bfile
      real(r8) :: minAngle, maxAngle
      type(L1BData_T) :: l1bField
      integer :: nomafs
      integer :: flag
      character(len=NameLen) :: l1bitemname
      real(r8) :: last, first
      real(r8) :: delta
      real(r8) :: incline
      real(r8), parameter :: secondsInDay = 24*60*60
      ! Note this next one is ONLY NEEDED for the case where we have only
      ! one MAF in the chunk
      real(r8), parameter :: ORBITALPERIOD = 98.8418*60.0
      real(r8), dimension(:), pointer :: mif1GeodAngle
      integer :: i ! Common counter for loops
      integer :: n, extra, left, right
      real(r8), dimension(:), pointer :: TMPANGLE ! A temporary array for the single case

      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      call ChooseOptimumLon0( filedatabase, chunk, 'GHz' )
      l1bItemName = AssembleL1BQtyName (instrumentModuleName // ".tpGeodAngle", &
                                        l1bfile%hdfVersion, .false.)
      call ReadL1BData (L1BFile, trim(l1bItemName), l1bField, noMafs, flag, &
                        firstMaf=chunk%firstMafIndex, lastMaf=chunk%lastMafIndex, &
                        dontPad=.true.)
      nullify (mif1GeodAngle)
      call Allocate_Test (mif1geodangle, nomafs, "mif1geodangle", modulename)
      mif1GeodAngle = l1bField%dpField(1,1,1:noMAFs)

      maxAngle = maxval (l1bField%dpField(1,:,nomafs))

      ! Now choose the geodetic angles for the hGrid
      mif1geodangle = l1bfield%dpField(1,1,1:nomafs)
      minAngle = minval(l1bField%dpField(1,:,1))
      first = origin + spacing * int ((minAngle-origin)/spacing)
      delta = first - minAngle ! it means first could be smaller
      if (delta > spacing/2) then
         first = first - spacing
      else if (delta < -spacing/2) then
         first = first + spacing
      end if

      ! Now work out the last point in a similar manner
      last = origin + spacing * int ((maxAngle-origin)/spacing)
      delta = last - maxAngle ! So +ve means last could be smaller
      if (delta > spacing/2) then
         last = last - spacing
      else if (delta < -spacing/2) then
         last = last + spacing
      end if

      call DeallocateL1BData ( l1bField )

      ! single is a non-existent concept for CFM
      ! Because the GMAO people, don't know of the existence of profiles
      if (.not. insetOverlaps) then
         if (first > minAngle) first = first - spacing
         if (last < maxAngle) last = last + spacing
      end if

      if ( chunk%firstMafIndex >= chunk%lastMafIndex ) last = first
      hgrid%noProfs = nint ( (last - first) / spacing ) + 1
      call CreateEmptyHGrid(hgrid)
      do i = 1, hgrid%noProfs
         hGrid%phi(1,i) = first + (i-1)*spacing
      end do

      ! Now fill the other geolocation information, first latitude
      ! Get orbital inclination
      l1bItemName = AssembleL1BQtyName ( "scOrbIncl", l1bfile%hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, &
          & l1bField, noMAFs, flag, &
          & firstMAF=chunk%firstMAFIndex, &
          & lastMAF=chunk%lastMAFIndex, &
          & dontPad=.true. )
      incline = sum(l1bField%dpField(1,1,:)) / nomafs
      call deallocateL1bData(l1bfield)

      hgrid%geodlat = rad2deg * asin (sin(deg2rad*hgrid%phi) * sin(deg2rad*incline))

      call EmpiricalLongitude(hgrid%phi(1,:), hgrid%lon(1,:))

      l1bitemname = AssembleL1BQtyName("MAFStartTimeTAI", l1bfile%hdfversion, .false.)
      call ReadL1BData(l1bfile, l1bItemName, l1bfield, nomafs, &
                       flag, firstMaf=chunk%firstMafIndex, &
                       lastMaf=chunk%lastMafIndex, dontpad=.true.)
      if (chunk%firstMafIndex /= chunk%lastMafIndex) then
         call InterpolateValues(mif1GeodAngle, l1bfield%dpField(1,1,:), &
                                hgrid%phi(1,:), hgrid%time(1,:), &
                                method='Linear', extrapolate='Allow')
      else
         ! Case where only single MAF per chunk, treat it specially
         hgrid%time = l1bfield%dpfield(1,1,1) + &
                      (OrbitalPeriod / 360.0) * &
                      (hgrid%phi - mif1GeodAngle(1))
      end if
      call DeallocateL1BData(l1bfield)
      ! First get fractional day, note this neglects leap seconds.
      ! Perhaps fix this later !???????? NJL. We do have access to the
      ! UTC ascii time field, perhaps we could use that?
      hgrid%solartime = modulo (hgrid%time, secondsInday) / secondsInDay
      ! Now correct for longitude and convert to hours
      hgrid%solarTime = 24.0 * (hgrid%solarTime + hGrid%lon/360.0)
      hgrid%solarTime = modulo (hGrid%solarTime, 24.0_r8)

      ! Solar zenith
      l1bItemName = AssembleL1BQtyName(instrumentModuleName//".tpSolarZenith", &
                                       l1bfile%hdfVersion, .false. )
      call ReadL1BData (l1bfile, l1bitemname, l1bfield, nomafs, &
                        flag, firstMaf=chunk%firstMafIndex, &
                        lastMaf=chunk%lastMafIndex, dontpad=.true.)
      call InterpolateValues (mif1GeodAngle, l1bField%dpField(1,1,:), &
                              hgrid%phi(1,:), hgrid%solarZenith(1,:), &
                              method='Linear', extrapolate='Allow')
      call DeallocateL1BData(l1bfield)

      ! Line of sight angle
      ! This we'll have to do with straight interpolation
      l1bitemname = AssembleL1BQtyName(instrumentModuleName//".tpLosAngle", &
                                       l1bfile%hdfVersion, .false.)
      call ReadL1BData(l1bfile, l1bitemname, l1bfield, nomafs, flag, &
                       firstMaf=chunk%firstMafIndex, &
                       lastMaf=chunk%lastMafIndex, dontPad=.true.)
      call InterpolateValues (mif1geodangle, l1bfield%dpField(1,1,:), &
                              hgrid%phi(1,:), hGrid%losAngle(1,:), &
                              method='Linear', extrapolate='Allow')
      call DeallocateL1BData(l1bfield)
      hgrid%losAngle = modulo(hgrid%losAngle, 360.0_r8)

      ! Done
      call Deallocate_test(mif1GeodAngle, 'mif1GeodAngle', modulename)

   end function

!--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
!---------------------------------------------------------------------------

end module

! $Log$
! Revision 1.12  2016/07/06 22:11:16  pwagner
! We seem to ned this change sometimes
!
! Revision 1.11  2011/12/15 23:47:49  honghanh
! Add documentation for CreateStandardHGrid
!
! Revision 1.10  2011/10/19 16:06:13  honghanh
! Add function to create hgrid with start maf and end maf.
!
! Revision 1.9  2010/06/29 17:02:47  honghanh
! Change the identifier 'fakeChunk' to 'chunk' because
! since it is created with ChunkDivide, it's as real as a chunk
! can get.
!
! Revision 1.8  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.7  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
