module CFM_Fill_m

   use Init_Tables_Module, only: L_LOSVEL, &
       L_L1BMIF_TAI, L_L1BMAFBASELINE, &
       L_ECRTOFOV, L_PTAN, L_ORBITINCLINATION, &
       L_RADIANCE, L_SCGEOCALT, L_SCVEL, &
       L_TNGTGEODALT, L_TNGTGEOCALT, L_SCVELECR, &
       L_TNGTECI, L_SCVELECI, L_SCECI
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error
   use MLSFiles, only: GetMLSFileByType, HDFVERSION_5
   use VectorsModule, only: VectorValue_T, M_LINALG, M_fill, MaskVectorQty
   use L1BData, only: GetL1BFile, ASSEMBLEL1BQTYNAME, L1BData_T, &
             DeallocateL1BData, ReadL1BData
   use MLSCommon, only: MLSFile_T, DEFAULTUNDEFINEDVALUE, r8
   use MLSSignals_m, only: GetSignalName, GetModuleName
   use BitStuff, only: NegativeIfBitPatternSet
   use Chunks_m, only: MLSChunk_T
   use MLSStrings, only: writeIntsToChars
   use FillUtils_1, only: fillerror, &
      Orig_FillVectorQuantityFromL1B => FillVectorQuantityFromL1B, &
      Orig_FillPhiTanWithRefraction => FillPhiTanWithRefraction

   implicit none

   public :: ExplicitFillVectorQuantity, FillVectorQuantityFromL1B
   public :: SpreadFillVectorQuantity
   public :: FillVectorQtyFromProfile

   private

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   contains

   ! Fill the quantity with given values.
   subroutine ExplicitFillVectorQuantity (quantity, values)
      ! The quantity to be filled. Only values are filled,
      ! other fields of VectorValue_T object won't be overwritten.
      type(VectorValue_T), intent(inout) :: quantity
      ! the amount of values provided must be equal to
      ! quantity%template%instanceLen * quantity%template%noInstances
      real(r8), dimension(:), intent(in) :: values

      integer :: noValues, numChans
      integer :: i,j,k
      integer :: surf, chan
      character(len=10) :: int1 = "          ", int2 = "          "

      noValues = size(values)

      if (noValues /= quantity%template%instanceLen * &
                quantity%template%noInstances) then
         call writeIntsToChars( &
            quantity%template%instanceLen * quantity%template%noInstances, int1)
         call writeIntsToChars(noValues, int2)
         call MLSMessage (MLSMSG_Error, moduleName, &
           "Incorrect amount of data, expect " // trim(int1) // ", has " // trim(int2))
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

   ! Fill all values in the quantity with the same number
   subroutine SpreadFillVectorQuantity (quantity, value)
      type(VectorValue_T), intent(inout) :: quantity
      real(r8), intent(in) :: value

      integer :: numChans
      integer :: k, i, j
      integer :: surf, chan

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
               quantity%values(j,i) = value
            end do
         end do
      end do
   end subroutine

   ! Fill the quantity with data from L1B files
   subroutine FillVectorQuantityFromL1B (quantity, chunk, filedatabase, &
      isPrecision, suffix, precisionQuantity, BOMask)
      ! The quantity to be filled
      type(VectorValue_T), intent(inout) :: quantity
      ! fake chunk object carrying needed information
      type(MLSChunk_T), intent(in) :: chunk
      ! the filedatabase containing L1B files
      type(MLSFile_T), dimension(:), pointer :: filedatabase
      ! is this a precision quantity
      logical, intent(in) :: isPrecision
      ! for now, don't use the optional arguments
      integer, intent(in), optional :: suffix
      type(VectorValue_T), intent(in), optional :: precisionQuantity
      integer, intent(in), optional :: BOMask

      fillError = 0
      call Orig_FillVectorQuantityFromL1B(0, quantity, chunk, filedatabase, &
         isPrecision, suffix, precisionQuantity, BOMask)
      if (fillError /= 0) then
         call MLSMessage (MLSMSG_Error, moduleName, "Can't Fill from L1B")
      end if
   end subroutine

   ! Don't use
   subroutine FillPhiTanWithRefraction (quantity, h2o, orbIncl, ptan, &
                                        refGPH, temperature)
      type(VectorValue_T), intent(inout) :: quantity
      type(VectorValue_T), intent(in) :: h2o
      type(VectorValue_T), intent(in) :: orbincl
      type(VectorValue_T), intent(in) :: ptan
      type(VectorValue_T), intent(in) :: refGPH
      type(VectorValue_T), intent(in) :: temperature

      call Orig_FillPhiTanWithRefraction(0, quantity, h2o, orbIncl, &
                                         ptan, refGPH, temperature)
   end subroutine

   ! Don't know how to describe this
   subroutine FillVectorQtyFromProfile (quantity, dontMask, heights, &
                                        values, value_unit, ptan, logSpace)
      use VectorsModule, only: ValidateVectorQuantity
      use VGridsDatabase, only: GETUNITFORVERTICALCOORDINATE
      use Init_Tables_Module, only: l_pressure, phyq_dimensionless, phyq_zeta
      use MLSNumerics, only: InterpolateValues, Hunt
      use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

      type(VectorValue_T), intent(inout) :: quantity
      logical, intent(in) :: dontMask
      real(r8), dimension(:), intent(in) :: heights
      real(r8), dimension(:), intent(in) :: values
      integer :: value_unit
      type(VectorValue_T), intent(in), optional :: ptan
      logical, intent(in), optional :: logSpace

      integer :: heightUnit, noUnique
      logical :: fail
      integer :: i, j, s, c, noPoints
      logical :: mylogSpace
      logical, dimension(:), pointer :: duplicated
      logical :: LOCALOUTHEIGHTS ! Set if out heights is our own variable
      real (r8), dimension(:), pointer :: OUTHEIGHTS ! Heights for output
      integer, dimension(:), pointer :: ININDS ! Indices
      real (r8), dimension(:), pointer :: OUTVALUES ! Single profile for output
      real(r8), dimension(size(heights)) :: myheights
      real(r8), dimension(size(values)) :: myvalues

      if (.not. ValidateVectorQuantity(quantity, coherent=.true.) &
          .and. .not. present(ptan)) &
          call MLSMessage(MLSMSG_Error, ModuleName, &
          'The quantity is not amenable to a profile fill unless you supply ptan')

      heightUnit = GetUnitForVerticalCoordinate(quantity%template%verticalCoordinate)
      if (present(ptan)) heightUnit = phyq_zeta

      myLogSpace = quantity%template%logBasis
      if (present(logSpace)) myLogSpace = logSpace

      if (value_unit /= phyq_dimensionless &
          .and. value_unit /= quantity%template%unit) &
          call MLSMessage(MLSMSG_Error, moduleName, "Bad unit for profile fill")

      if (size(heights) /= size(values)) &
         call MLSMessage(MLSMSG_Error, moduleName, &
         "Mismatch length between heights and values array")

      noPoints = size(heights)
      nullify(duplicated, outheights, outvalues, ininds)
      call Allocate_test ( duplicated, noPoints, 'duplicated', ModuleName )
      call Allocate_test (outValues, quantity%template%noSurfs, 'outValues', ModuleName)

      if (heightUnit == phyq_zeta) then
         myheights = -log10(heights)
      else
         myheights = heights
      end if

      if (myLogSpace .and. any(values <= 0.0)) &
         call MLSMessage(MLSMSG_Error, moduleName, &
         'Non-positive input data in log profile fill')

      if (myLogSpace) then
         myValues = log(values)
      else
         myValues = values
      end if

      ! Get the appropriate height coordinate for output, for pressure take log.
      if ( present(ptan) ) then
        localOutHeights = .false.
        outHeights => ptan%values(:,1)
      elseif ( quantity%template%verticalCoordinate == l_pressure ) then
        localOutHeights = .true.
        call Allocate_test (outHeights, quantity%template%noSurfs, &
                            'outHeights', ModuleName )
        outHeights = -log10 ( quantity%template%surfs(:,1) )
      else
        localOutHeights = .false.
        outHeights => quantity%template%surfs(:,1)
      end if

      ! Now, if the quantity is coherent, let's assume the user wanted the
      ! 'nearest' values
      if ( quantity%template%coherent .or. present(ptan) ) then
        nullify ( inInds )
        call allocate_test ( inInds, noPoints, 'inInds', ModuleName )
        call hunt ( outHeights, myHeights, inInds, &
          & nearest=.true., allowTopValue=.true., fail=fail )
        if ( fail ) then
          call MLSMessage ( MLSMSG_Error, moduleName, 'Problem in Hunt' )
        end if
        duplicated = .false.
        do i = 1, noPoints - 1
          do j = i + 1, noPoints
            if ( inInds(i) == inInds(j) ) then
              duplicated ( j ) = .true.
            end if
          end do
        end do
        noUnique = count ( .not. duplicated )
        inInds(1:noUnique) = pack ( inInds, .not. duplicated )
        myHeights(1:noUnique) = outHeights ( inInds(1:noUnique) )
        myValues(1:noUnique) = pack ( myValues, .not. duplicated )
        call deallocate_test ( inInds, 'inInds', ModuleName )
      end if

      ! Now do the interpolation for the first instance
      call InterpolateValues ( myHeights(1:noUnique), myValues(1:noUnique), outHeights, &
        & outValues, 'Linear', extrapolate='Constant' )

      if (myLogSpace) outvalues = exp(outvalues)

      do i = 1, quantity%template%noInstances
         j = 1
         do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
               if (associated(quantity%mask) .and. .not. dontMask) then
                  if (iand(ichar(quantity%mask(j,i)), m_fill) == 0) &
                     quantity%values(j,i) = outValues(s)
               else
                  quantity%values(j,i) = outValues(s)
               end if
               j = j+1
            end do
         end do
      end do

      ! Finish off
      if ( localOutHeights ) call Deallocate_test ( outHeights, &
         & 'outHeights', ModuleName )
      call Deallocate_test ( duplicated, 'duplicated', ModuleName )
      call Deallocate_test ( outValues, 'outValues', ModuleName )

   end subroutine

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
