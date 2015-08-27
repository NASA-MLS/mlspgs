! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_Fill_m

    use INIT_TABLES_MODULE, only: L_LOSVEL, &
        L_L1BMIF_TAI, L_L1BMAFBASELINE, L_NONE, &
        L_ECRTOFOV, L_PTAN, L_ORBITINCLINATION, &
        L_RADIANCE, L_SCGEOCALT, &
        L_TNGTGEODALT, L_TNGTGEOCALT, L_SCVELECR, &
        L_TNGTECI, L_SCVELECI, L_SCECI
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSFILES, only: GETMLSFILEBYTYPE, HDFVERSION_5
    use VECTORSMODULE, only: VECTORVALUE_T, M_LINALG, M_FILL, MASKVECTORQTY, DUMP
    use L1BDATA, only: GETL1BFILE, ASSEMBLEL1BQTYNAME, L1BDATA_T, &
              DEALLOCATEL1BDATA, READL1BDATA
    use MLSCOMMON, only: MLSFILE_T, DEFAULTUNDEFINEDVALUE
    use MLSKINDS, only: R8
    use MLSSIGNALS_M, only: GETSIGNALNAME, GETMODULENAME
    use BITSTUFF, only: NEGATIVEIFBITPATTERNSET
    use CHUNKS_M, only: MLSCHUNK_T
    use MLSSTRINGS, only: WRITEINTSTOCHARS
    use FILLUTILS_1, only: FILLERROR, FROML1B, &
                           PHITANWITHREFRACTION
    use STRING_TABLE, only: CREATE_STRING

    implicit none

    public :: EXPLICITFILLVECTORQUANTITY, FILLVECTORQUANTITYFROML1B
    public :: SPREADFILLVECTORQUANTITY, FILLPTANQUANTITY
    public :: FILLVECTORQTYFROMPROFILE, FILLPHITANQUANTITY
    public :: APPLYBASELINE

    private

!---------------------------- RCS Ident Info -------------------------------
    character(len=*), private, parameter :: ModuleName= &
        "$RCSfile$"
!---------------------------------------------------------------------------

    interface FillVectorQuantityFromL1B
        module procedure FillVectorQuantityFromL1B_maf
        module procedure FillVectorQuantityFromL1B_chunk
    end interface

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

    ! Fill the quantity with data from L1B files with maf range
    ! specified by firstL1Maf and lastL1Maf
    subroutine FillVectorQuantityFromL1B_maf (quantity, firstL1Maf, lastL1Maf, &
    filedatabase, isPrecision, suffix, precisionQuantity, BOMask)
        ! The quantity to be filled
        type(VectorValue_T), intent(inout) :: quantity
        ! The maf range of the data to be read from L1B files
        integer, intent(in) :: firstL1Maf, lastL1Maf
        ! the filedatabase containing L1B files
        type(MLSFile_T), dimension(:), pointer :: filedatabase
        ! is this a precision quantity
        logical, intent(in) :: isPrecision
        character(len=*), intent(in), optional :: suffix
        type(VectorValue_T), intent(in), optional :: precisionQuantity
        ! If isPrecision is .true. and BOMask is not 0, then
        ! bright object status is read from L1BOA file
        integer, intent(in), optional :: BOMask

        type (MLSChunk_T) :: Chunk

        chunk%firstMafIndex = firstL1Maf
        chunk%lastMafIndex = lastL1Maf

        call FillVectorQuantityFromL1B_chunk (quantity, chunk, filedatabase, &
        isPrecision, suffix, precisionQuantity, BOMask)
    end subroutine

    ! Fill the quantity with data from L1B files
    subroutine FillVectorQuantityFromL1B_chunk (quantity, chunk, filedatabase, &
    isPrecision, suffix, precisionQuantity, BOMask)
        ! The quantity to be filled
        type(VectorValue_T), intent(inout) :: quantity
        ! fake chunk object carrying needed information
        type(MLSChunk_T), intent(in) :: chunk
        ! the filedatabase containing L1B files
        type(MLSFile_T), dimension(:), pointer :: filedatabase
        ! is this a precision quantity
        logical, intent(in) :: isPrecision
        character(len=*), intent(in), optional :: suffix
        type(VectorValue_T), intent(in), optional :: precisionQuantity
        ! If isPrecision is .true. and BOMask is not 0, then
        ! bright object status is read from L1BOA file
        integer, intent(in), optional :: BOMask
        integer :: geoLocation
        integer :: i

        fillError = 0
        geoLocation = l_none
        i = 0
        if (present (suffix)) then
            i = create_string(suffix)
        end if

        call FromL1B(0, quantity, chunk, filedatabase, &
                     isPrecision, i, geoLocation, precisionQuantity, BOMask)
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

        call PhiTanWithRefraction( 0, quantity, h2o, orbIncl, &
                                  ptan, refGPH, temperature, .false. )
    end subroutine

    ! Don't know how to describe this
    subroutine FillVectorQtyFromProfile (quantity, dontMask, heights, &
                                         values, value_unit, ptan, logSpace)
        use VectorsModule, only: ValidateVectorQuantity
        use VGridsDatabase, only: GETUNITFORVERTICALCOORDINATE
        use Init_Tables_Module, only: l_pressure, phyq_dimensionless, phyq_zeta
        use MLSNumerics, only: InterpolateValues, Hunt
        use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

        type(VectorValue_T), intent(inout), target :: quantity
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

    ! Fill phitan quantity with phi values from the quantity's template
    subroutine FillPhitanQuantity (quantity)
        ! phitan quantity
        type(VectorValue_T), intent(inout) :: quantity

        quantity%values = quantity%template%phi
    end subroutine

    ! Compute and fill ptan quantity from other quantities
    subroutine FillPtanQuantity (ptan, temperature, refGPH, h2o, orbitInclination, &
                                 phitan, geocentricAltitude)
        use ScanModelModule, only: Get2DHydrostaticTangentPressure
        use Init_tables_module, only: phyq_angle

        type(VectorValue_T), intent(inout) :: ptan
        ! Input quantities
        type(VectorValue_T), intent(in) :: temperature, refGPH
        type(VectorValue_T), intent(in) :: h2o, orbitInclination, phitan
        type(VectorValue_T), intent(in) :: geocentricAltitude

        call Get2DHydrostaticTangentPressure(ptan, temperature, refGPH, &
        & h2o, orbitInclination, phitan, geocentricAltitude, 4, &
        & (/ 0.0_r8, 0.0_r8 /), phyq_angle)
    end subroutine

    ! Applies baseline quantity to a radiance quantity or applies noise
    ! to a precision quantity of a radiance quantity.
    subroutine ApplyBaseline (quantity, baselineQuantity, quadrature, dontmask)
        use FillUtils_1, only: Orig_ApplyBaseline => ApplyBaseline

        ! Radiance quantity to modify
        type(VectorValue_T), intent(inout) :: quantity
        ! L1B MAF baseline to use
        type(VectorValue_T), intent(in) :: baselineQuantity
        logical, intent(in) :: quadrature
        logical, intent(in) :: dontmask

        call Orig_ApplyBaseline ( 0, quantity, baselineQuantity, quadrature, &
          & dontmask, .false. )
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

! $Log$
! Revision 1.17  2015/08/05 20:21:29  pwagner
! Modified to compile properly with v4
!
! Revision 1.16  2013/07/26 16:48:15  pwagner
! Consistent with removal of SCVEL
!
! Revision 1.15  2013/07/10 17:51:24  pwagner
! Changed to be consistent with FromL1B api
!
! Revision 1.14  2012/10/15 17:11:41  pwagner
! Adapted to new api in FillUtils
!
! Revision 1.13  2011/12/15 18:27:44  honghanh
! Documentation and code clean up, including removing unused and broken
! subroutines.
!
! Revision 1.12  2011/11/03 14:39:57  honghanh
! Add fill subroutine that use maf instead of chunk
!
! Revision 1.11  2010/07/08 21:39:16  honghanh
! Add ApplyBaseline to cfm_fill_m
!
! Revision 1.9  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
