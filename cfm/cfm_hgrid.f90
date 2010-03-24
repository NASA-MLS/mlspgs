module CFM_HGrid
   use HGridsDatabase, only: HGrid_T, CreateEmptyHGrid, TrimHGrid, &
                             FindClosestMatch, DestroyHGridContents, Dump
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

   public :: CreateRegularHGrid, DestroyHGridContents, Dump
   public :: HGrid_T

   private
   character(len=20), parameter :: moduleName = "CFM_HGrid"

   contains
   type(HGrid_T) function CreateRegularHGrid (instrumentModuleName, origin, &
                  spacing, insetOverLaps, single, filedatabase, &
                  fakeChunk) result(hgrid)

      real(r8), intent(in) :: origin
      real(r8), intent(in) :: spacing
      type(MLSFile_T), dimension(:), pointer :: filedatabase
      character(len=*), intent(in) :: instrumentModuleName
      logical, intent(in) :: insetOverlaps
      logical, intent(in) :: single
      type(MLSChunk_T) :: fakeChunk

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

      ! Executables
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      call ChooseOptimumLon0(filedatabase, fakeChunk)
      l1bItemName = AssembleL1BQtyName (instrumentModuleName // ".tpGeodAngle", &
                                        l1bfile%hdfVersion, .false.)
      call ReadL1BData (L1BFile, trim(l1bItemName), l1bField, noMafs, flag, &
                        firstMaf=fakeChunk%firstMafIndex, lastMaf=fakeChunk%lastMafIndex, &
                        dontPad=.true.)
      nullify (mif1GeodAngle)
      call Allocate_Test (mif1geodangle, nomafs, "mif1geodangle", modulename)
      mif1GeodAngle = l1bField%dpField(1,1,1:noMAFs)

      maxAngle = maxval (l1bField%dpField(1,:,nomafs))

      ! Now choose the geodetic angles for the hGrid
      if (.not. single) then
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
      else
         ! The 'single' option is typically used for running single profile retrievals
         ! in a debug mode.  In order to ensure we choose the same profile for each MAF
         ! that a phiWindow=0 forward model would select, we have to do some extra work.

         ! First check this is a sane request
         if ( noMAFs /= 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
           & 'Single hGrid option set but more than one MAF in chunk' )

         ! Now construct a temporary HGrid that spans an orbit either side of this MAF
         ! Overkill I know but not harmful
         first = origin + spacing * int ( (minAngle-origin)/spacing ) - 360.0_r8
         last = first + 720.0_r8
         n = (last - first) / spacing
         nullify(tmpAngle)
         call Allocate_test ( tmpAngle, n, 'tmpAngle', ModuleName )
         do i = 1,n
            tmpAngle(i) = first + (i-1) * spacing
         end do
         i = FindClosestMatch ( tmpAngle, l1bField%dpField(1,:,:), 1 )
         first = tmpAngle(i)
         last = first
         call Deallocate_test ( tmpAngle, 'tmpAngle', ModuleName )
      end if

      call DeallocateL1BData ( l1bField )

      if (.not. insetOverlaps .and. .not. single) then
         if (first > minAngle) first = first - spacing
         if (last < maxAngle) last = last + spacing
      end if

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
          & firstMAF=fakeChunk%firstMAFIndex, &
          & lastMAF=fakeChunk%lastMAFIndex, &
          & dontPad=.true. )
      incline = sum(l1bField%dpField(1,1,:)) / nomafs
      call deallocateL1bData(l1bfield)

      hgrid%geodlat = rad2deg * asin (sin(deg2rad*hgrid%phi) * sin(deg2rad*incline))

      call EmpiricalLongitude(hgrid%phi(1,:), hgrid%lon(1,:))

      l1bitemname = AssembleL1BQtyName("MAFStartTimeTAI", l1bfile%hdfversion, .false.)
      call ReadL1BData(l1bfile, l1bItemName, l1bfield, nomafs, &
                       flag, firstMaf=fakeChunk%firstMafIndex, &
                       lastMaf=fakeChunk%lastMafIndex, dontpad=.true.)
      if (fakeChunk%firstMafIndex /= fakeChunk%lastMafIndex) then
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
                        flag, firstMaf=fakeChunk%firstMafIndex, &
                        lastMaf=fakeChunk%lastMafIndex, dontpad=.true.)
      call InterpolateValues (mif1GeodAngle, l1bField%dpField(1,1,:), &
                              hgrid%phi(1,:), hgrid%solarZenith(1,:), &
                              method='Linear', extrapolate='Allow')
      call DeallocateL1BData(l1bfield)

      ! Line of sight angle
      ! This we'll have to do with straight interpolation
      l1bitemname = AssembleL1BQtyName(instrumentModuleName//".tpLosAngle", &
                                       l1bfile%hdfVersion, .false.)
      call ReadL1BData(l1bfile, l1bitemname, l1bfield, nomafs, flag, &
                       firstMaf=fakeChunk%firstMafIndex, &
                       lastMaf=fakeChunk%lastMafIndex, dontPad=.true.)
      call InterpolateValues (mif1geodangle, l1bfield%dpField(1,1,:), &
                              hgrid%phi(1,:), hGrid%losAngle(1,:), &
                              method='Linear', extrapolate='Allow')
      call DeallocateL1BData(l1bfield)
      hgrid%losAngle = modulo(hgrid%losAngle, 360.0_r8)

      if (single) then
         ! Delete all but the 'center' profile
         extra = hGrid%noProfs - 1
         left = extra / 2
         right = extra - left
         if (left > 0) call TrimHGrid(hGrid, -1, left)
         if (right > 0) call TrimHGrid(hGrid, -1, right)
      end if !We know that upper/lower overlap is guaranteed to be 0

      ! Done
      call Deallocate_test(mif1GeodAngle, 'mif1GeodAngle', modulename)

   end function

end module
