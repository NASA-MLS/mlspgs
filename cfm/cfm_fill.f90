module CFM_Fill
 
   use Init_Tables_Module, only: L_LOSVEL, &
       L_L1BMIF_TAI, L_L1BMAFBASELINE, &
       L_ECRTOFOV, L_PTAN, L_ORBITINCLINATION, &
       L_RADIANCE, L_SCGEOCALT, L_SCVEL, &
       L_TNGTGEODALT, L_TNGTGEOCALT, L_SCVELECR, &
       L_TNGTECI, L_SCVELECI, L_SCECI
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error
   use MLSFiles, only: GetMLSFileByType, HDFVERSION_5
   use VectorsModule, only: VectorValue_T, M_LINALG, MaskVectorQty
   use L1BData, only: GetL1BFile, ASSEMBLEL1BQTYNAME, L1BData_T, &
             DeallocateL1BData, ReadL1BData
   use MLSCommon, only: MLSFile_T, DEFAULTUNDEFINEDVALUE, r8
   use MLSSignals_m, only: GetSignalName, GetModuleName
   use BitStuff, only: NegativeIfBitPatternSet
   use Chunks_m, only: MLSChunk_T
  
   implicit none

   public :: ExplicitFillVectorQuantity

   private
   character(len=20) :: moduleName="CFM_Fill"

   contains
   subroutine FillVectorQuantityFromL1B (quantity, chunk, filedatabase)
      type(VectorValue_T), intent(inout) :: quantity
      type(MLSChunk_T), intent(in) :: chunk
      type(MLSFile_T), dimension(:), pointer :: filedatabase

      type(MLSFile_T), pointer :: l1bFile
      integer :: hdf_version
      character(len=132) :: namestring
      type(L1BData_T) :: data, BO_stat
      integer :: flag, noMafs, myBOMask, channel, maxMifs
      integer :: row, column

      myBOMask = 0
      l1bfile => GetMLSFileByType (filedatabase, content='l1boa')
      hdf_version = l1bfile%hdfVersion

      select case (quantity%template%quantityType)
      case (l_ECRtoFOV)
         call GetModuleName (quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('ECRtoFOV', hdf_version, &
              .true., trim(nameString))
      case (l_l1bMAFBaseline)
         call GetSignalName (quantity%template%signal, nameString, &
            sideband=quantity%template%sideband, noChannels=.true.)
         nameString = AssembleL1BQtyName(nameString, hdf_version, .false.)
      case (l_l1bMIF_TAI)
         if (hdf_version == HDFVERSION_5) then
            call GetModuleName (quantity%template%instrumentModule, nameString)
            nameString = AssembleL1BQtyName('MIF_TAI', hdf_version, .false., &
                 trim(nameString))
         else 
            nameString = 'MIF_TAI'
         end if
      case (l_losVel)
         call GetModuleName(quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('LosVel', hdf_version, .true., &
              trim(nameString))
      case (l_orbitInclination)
         nameString = AssembleL1BQtyName('OrbIncl', hdf_version, .false., &
              'sc')
      case (l_ptan)
         call GetModuleName (quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('ptan', hdf_version, .false., &
             trim(nameString))
      case (l_radiance)
         call GetSignalName(quantity%template%signal, nameString, &
             sideband=quantity%template%sideband, noChannels=.true.)
         nameString = AssembleL1BQtyName(nameString, hdf_version, .false.)
         l1bfile => GetL1BFile(filedatabase, nameString)
      case (l_scECI)
         nameString = AssembleL1BQtyName('ECI', hdf_version, .false., 'sc')
      case (l_scGeocAlt)
         nameString = AssembleL1BQtyName('GeocAlt', hdf_version, .false., &
            'sc')
      case (l_scVel)
         nameString = AssembleL1BQtyName('Vel', hdf_version, .false., 'sc')
      case (l_scVelECI)
         nameString = AssembleL1BQtyName('VelECI', hdf_version, .false., 'sc')
      case (l_scVelECR)
         nameString = AssembleL1BQtyName('VelECR', hdf_version, .false., 'sc')
      case (l_tngtECI)
         call GetModuleName(quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('ECI', hdf_version, .true., &
            trim(nameString))
      case (l_tngtGeocAlt)
         call GetModuleName(quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('GeocAlt', hdf_version, .true., &
             trim(nameString))
      case (l_tngtGeodAlt)
         call GetModuleName(quantity%template%instrumentModule, nameString)
         nameString = AssembleL1BQtyName('GeodAlt', hdf_version, .true., &
            trim(nameString))
      case default
         print *, "Unknown quantity type"
         call MLSMessage(MLSMSG_Error, ModuleName, 'Unknown quantity type')
      end select

      if ( associated(L1BFile) .or. ( quantity%template%quantityType /= l_radiance .and. &
        & quantity%template%quantityType /= l_L1BMAFBaseline ) ) then
        L1BFile => GetL1bFile(filedatabase, namestring)
        call ReadL1BData ( L1BFile, nameString, data, noMAFs, flag, &
          & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & NeverFail= .false., &
          & dontPad=.false. )
        ! If it didn't exist in the not-a-radiance case, then we'll fail here.
        if (flag /= 0) then
           call MLSMessage(MLSMSG_Error, ModuleName, 'Error in reading L1B data')
        end if

       ! do channel = 1, quantity%template%noChans
       !   data%dpField(channel,:,:) = &
       !     & NegativeIfBitPatternSet( data%dpField(channel,:,:), &
       !     & BO_stat%intField(1, 1:maxMIFs, 1:noMAFs), myBOMask )
       ! enddo
       quantity%values = RESHAPE(data%dpField, &
          & (/ quantity%template%instanceLen, quantity%template%noInstances /) )
      else
       ! This is the case where it's a radiance we're after and it's missing
       quantity%values = DEFAULTUNDEFINEDVALUE ! -1.0
       do column=1, size(quantity%values(1,:))
          do row=1, size(quantity%values(:,1))
            call MaskVectorQty ( quantity, row, column, M_LinAlg )
          end do
       end do
      end if

   end subroutine

   subroutine ExplicitFillVectorQuantity (quantity, values)
      type(VectorValue_T), intent(inout) :: quantity
      real(r8), dimension(:), intent(in) :: values

      integer :: noValues, numChans
      integer :: i,j,k
      integer :: surf, chan

      noValues = size(values)

      if (noValues /= quantity%template%instanceLen * &
                quantity%template%noInstances) then
         ! need to find out how to convert integer to char to write
         ! better error message here
         print *, "not the right amount of data, expect ", &
                  quantity%template%instanceLen * quantity%template%noInstances, &
                  ", has ", noValues
         call MLSMessage (MLSMSG_Error, moduleName, &
           "Incorrect amount of data")
      end if

      ! need checking on the value and their units?

      numChans = quantity%template%instanceLen / quantity%template%noSurfs
      if (numChans /= quantity%template%noChans) then
         call MLSMessage (MLSMSG_Error, moduleName, &
           "Inconsistent template instance length")
      end if

      ! loop thru the quantity
      k = 0
      do i = 1, quantity%template%noInstances
         j = 0
         do surf = 1, quantity%template%noSurfs
            do chan = 1, numChans
               j = j + 1
               k = k + 1
               quantity%values(j,i) = values (mod (k-1, noValues) + 1)
            end do
         end do
      end do
   end subroutine
end module
