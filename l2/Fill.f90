! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
  !=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Expr_M, only: EXPR
  use GriddedData, only: GriddedData_T
  use INIT_TABLES_MODULE, only: F_GEOCALTITUDEQUANTITY, F_EXPLICITVALUES, &
    & F_H2OQUANTITY, F_MAXITERATIONS, F_METHOD, F_QUANTITY, F_REFGPHQUANTITY, &
    & F_SCECI, F_SCVEL, F_SOURCE, F_SOURCEGRID, F_SOURCEL2AUX, F_SOURCEL2GP, &
    & F_SOURCEQUANTITY, F_SPREAD, F_TEMPERATUREQUANTITY, F_TNGTECI, FIELD_FIRST, &
    & FIELD_LAST, L_EXPLICIT, L_GPH, L_GRIDDED, L_HYDROSTATIC, L_L1B, L_L2GP, &
    & L_L2AUX, L_LOSVEL, L_PRESSURE, L_PTAN, L_RADIANCE, L_REFGPH, &
    & L_SCECI, L_SCGEOCALT, L_SCVEL, L_SPECIAL, L_TEMPERATURE, L_TNGTECI,&
    & L_TNGTGEODALT, L_TNGTGEOCALT, L_TRUE, L_VECTOR, L_VMR, &
    & L_ZETA, S_TIME, S_VECTOR, S_FILL, S_SNOOP
  ! will be added
  use L1BData, only: DeallocateL1BData, FindL1BData, L1BData_T, ReadL1BData
  use L2GPData, only: L2GPData_T
  use L2AUXData, only: L2AUXData_T, L2AUXRank
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, NameLen, LineLen, MLSChunk_T, R8
  use MLSSignals_m, only: GetSignalName, GetModuleName
  use Molecules, only: L_H2O
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use string_table, only: get_string
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, &
    & SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_DOT, N_SET_ONE
  use VectorsModule, only: AddVectorToDatabase, CreateVector, Dump, &
    & GetVectorQtyByTemplateIndex, ValidateVectorQuantity, Vector_T, &
    & VectorTemplate_T, VectorValue_T
  use ScanModelModule, only: GetBasisGPH, GetHydrostaticTangentPressure, OMEGA
  use Intrinsic, only: L_CHANNEL, L_INTERMEDIATEFREQUENCY, L_USBFREQUENCY,&
    & L_LSBFREQUENCY, L_MIF, L_MAF, PHYQ_Dimensionless, PHYQ_Invalid
  use SnoopMLSL2, only: SNOOP

  implicit none
  private
  public :: MLSL2Fill

  ! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

  ! Error codes for "announce_error"  
  integer, parameter :: WRONG_NUMBER = 1     ! of fields of a VECTOR command
  integer, parameter :: unknownQuantityName = WRONG_NUMBER+1
  integer, parameter :: source_not_in_db = unknownQuantityName+1
  integer, parameter :: zeroProfilesFound = source_not_in_db+1
  integer, parameter :: zeroGeodSpan = zeroProfilesFound+1
  integer, parameter :: vectorWontMatchL2GP = zeroGeodSpan+1
  integer, parameter :: cantFillFromL2AUX = vectorWontMatchL2GP+1
  integer, parameter :: vectorWontMatchPDef = cantFillFromL2AUX+1
  integer, parameter :: cantFillFromL1B = vectorWontMatchPDef+1
  ! more Error codes relating to FillVector
  integer, parameter :: numInstancesisZero = cantFillFromL1B+1
  integer, parameter :: numSurfsisZero = numInstancesisZero+1
  integer, parameter :: numChansisZero = numSurfsisZero+1
  integer, parameter :: objIsFullRank3 = numChansisZero+1
  integer, parameter :: otherErrorInFillVector = objIsFullRank3+1
  integer, parameter :: noSourceGridGiven= otherErrorInFillVector+1
  integer, parameter :: noSourceL2GPGiven= noSourceGridGiven+1
  integer, parameter :: noSourceL2AUXGiven= noSourceL2GPGiven+1
  integer, parameter :: noExplicitValuesGiven= noSourceL2AUXGiven+1
  integer, parameter :: noSourceQuantityGiven= noExplicitValuesGiven+1
  integer, parameter :: invalidExplicitFill= noSourceQuantityGiven+1
  integer, parameter :: badUnitsForExplicit= invalidExplicitFill+1

  ! Error codes resulting from squeeze
  integer, parameter :: n1_is_zero = badUnitsForExplicit+1
  integer, parameter :: n2_is_zero = n1_is_zero+1
  integer, parameter :: n3_is_zero = n2_is_zero+1
  integer, parameter :: m1_too_small = n3_is_zero+1
  integer, parameter :: m2_too_small = m1_too_small+1
  integer, parameter :: not_permutation = m2_too_small+1
  integer, parameter :: allocation_err = not_permutation+1
  integer, parameter :: deallocation_err = allocation_err+1

  ! miscellaneous
  integer, parameter :: miscellaneous_err = deallocation_err+1
  integer, parameter :: errorReadingL1B = miscellaneous_err+1
  integer, parameter :: needTempREFGPH = errorReadingL1B+1
  integer, parameter :: needH2O = needTempRefGPH+1
  integer, parameter :: needGeocAltitude = needH2O+1
  integer, parameter :: badGeocAltitudeQuantity = needGeocAltitude+1
  integer, parameter :: badTemperatureQuantity = badGeocAltitudeQuantity+1
  integer, parameter :: badREFGPHQuantity = badTemperatureQuantity+1
  integer, parameter :: nonConformingHydrostatic = badREFGPHQuantity+1
  integer, parameter :: badUnitsForMaxIterations = nonConformingHydrostatic+1
  integer, parameter :: noSpecialFill = badUnitsForMaxIterations + 1
  integer, parameter :: badlosVelFill = noSpecialFill + 1

  !  integer, parameter :: s_Fill = 0   ! to be replaced by entry in init_tables_module
  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
    "$id: fill.f90,v 1.1 2000/01/21 21:04:06 livesey Exp $"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module performs the Fill operation in the Level 2 software.  
  ! This takes a vector template, and creates and fills an appropriate vector

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( root, l1bInfo, griddedData, vectorTemplates, vectors, &
    & qtyTemplates , L2GPDatabase, L2AUXDatabase, chunks, chunkNo )

    ! This is the main routine for the module.  It parses the relevant lines
    ! of the l2cf and works out what to do.

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (L1BInfo_T), intent(in) :: l1bInfo
    type (GriddedData_T), dimension(:), pointer :: griddedData
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (Vector_T), dimension(:), pointer :: vectors
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
    type (L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    integer, intent(in) :: chunkNo

    ! Local variables
    type (VectorValue_T), pointer :: QUANTITY ! Quantity to be filled

    type (VectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    type (VectorValue_T), pointer :: H2OQUANTITY
    type (VectorValue_T), pointer :: REFGPHQUANTITY
    type (VectorValue_T), pointer :: SCECIQUANTITY
    type (VectorValue_T), pointer :: SCVELQUANTITY
    type (VectorValue_T), pointer :: TEMPERATUREQUANTITY
    type (VectorValue_T), pointer :: TNGTECIQUANTITY

    type (Vector_T) :: newVector ! A vector we've created
    character (LEN=LineLen) ::              msr
    real :: T1, T2              ! for timing
    
    integer :: GEOCALTITUDEQUANTITYINDEX    ! In the source vector
    integer :: GEOCALTITUDEVECTORINDEX      ! In the vector database
    integer :: ERRORCODE                ! 0 unless error; returned by called routines
    integer :: FIELDINDEX               ! Entry in tree
    integer :: FILLMETHOD               ! How will we fill this quantity
    integer :: GSON                     ! Descendant of Son
    integer :: GRIDINDEX                ! Index of requested grid
    integer :: H2OQUANTITYINDEX         ! in the quantities database
    integer :: H2OVECTORINDEX           ! In the vector database
    integer :: I, J, K                  ! Loop indices for section, spec, expr
    integer :: ind                      ! Temoprary index
    integer :: KEY                      ! Definitely n_named
    integer :: L2AUXINDEX               ! Index into L2AUXDatabase
    integer :: L2GPINDEX                ! Index into L2GPDatabase
    integer :: L2INDEX                  ! Where source is among l2gp or l2aux database
    integer :: MAXITERATIONS            ! For hydrostatic fill
    integer :: PREVDEFDQT               !
    integer :: QUANTITYINDEX            ! Within the vector
    integer :: REFGPHQUANTITYINDEX      ! in the quantities database
    integer :: REFGPHVECTORINDEX        ! In the vector database
    integer :: SCECIVECTORINDEX         ! In the vector database
    integer :: SCECIQUANTITYINDEX       ! In the quantities database
    integer :: SCVELVECTORINDEX         ! In the vector database
    integer :: SCVELQUANTITYINDEX       ! In the quantities database
    integer :: SON                      ! Of root, an n_spec_args or a n_named
    integer :: SOURCEQUANTITYINDEX      ! in the quantities database
    integer :: SOURCEVECTORINDEX        ! In the vector database
    integer :: TEMPERATUREQUANTITYINDEX ! in the quantities database
    integer :: TEMPERATUREVECTORINDEX   ! In the vector database
    integer :: TEMPLATEINDEX            ! In the template database
    integer :: TNGTECIVECTORINDEX       ! In the vector database
    integer :: TNGTECIQUANTITYINDEX     ! In the quantities database
    integer, DIMENSION(2) :: UNITASARRAY ! From expr
    real(r8), DIMENSION(2) :: VALUEASARRAY ! From expr
    integer :: VALUESNODE               ! For the parser
    integer :: VECTORINDEX              ! In the vector database
    integer :: VECTORNAME               ! Name of vector to create
                                        !
    logical :: TIMING
    logical :: SPREAD           ! Do we spread values accross instances in explict

    logical, DIMENSION(field_first:field_last) :: got

    !    INTEGER :: OL2FileHandle

    ! Executable code
    timing = .false.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )

    ! Logical id of file(s) holding old L2GP data
    !    OL2FileHandle = mlspcf_ol2gp_start

    ! starting quantities number for *this* vector; what if we have more?
    !    qtiesStart = 1
    !   Calculate qtiesStart for the specific quantity below

    error = 0
    templateIndex = -1
    vectorIndex = -1
    spread=.FALSE.
    maxIterations = 4

    ! Loop over the lines in the configuration file

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        vectorName = sub_rosa(subtree(1,son))
      else
        key = son
        vectorName = 0
      end if
      got= .false.

      ! Node_id(key) is now n_spec_args.

      select case( get_spec_id(key) )
      case ( s_vector )
        if ( nsons(key) /= 2 ) call announce_error ( son, wrong_number )
        templateIndex = decoration(decoration(subtree(2,subtree(2,key))))

        ! Create the vector, and add it to the database.

        call decorate ( key, AddVectorToDatabase ( vectors, &
          & CreateVector ( vectorName, vectorTemplates(templateIndex), &
          & qtyTemplates ) ) )

        ! That's the end of the create operation

      case ( s_fill )
        ! Now we're on actual Fill instructions.
        ! Loop over the instructions to the Fill command

        do j=2,nsons(key)
          gson = subtree(j,key) ! The argument
          fieldIndex=decoration(subtree(1,gson))
          if (nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
          got(fieldIndex)=.TRUE.
          select case ( fieldIndex )
          case (f_quantity)   ! What quantity are we filling quantity=vector.quantity
            vectorIndex=decoration(decoration(subtree(1,gson)))
            quantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_method)   ! How are we going to fill it?
            fillMethod=decoration(gson)
          case (f_sourceQuantity)       ! When filling from a vector, what vector/quantity
            sourceVectorIndex=decoration(decoration(subtree(1,gson)))
            sourceQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_tngtECI)              ! For special fill of losVel
            tngtECIVectorIndex=decoration(decoration(subtree(1,gson)))
            tngtECIQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_scECI)                ! For special fill of losVel
            scECIVectorIndex=decoration(decoration(subtree(1,gson)))
            scECIQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_scVel)                ! For special fill of losVel
            scVelVectorIndex=decoration(decoration(subtree(1,gson)))
            scVelQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_sourceL2AUX)          ! Which L2AUXDatabase entry to use
            l2auxIndex=decoration(decoration(gson))
          case (f_sourceL2GP)           ! Which L2GPDatabase entry to use
            l2gpIndex=decoration(decoration(gson))
          case (f_temperatureQuantity) ! For hydrostatic
            temperatureVectorIndex=decoration(decoration(subtree(1,gson)))
            temperatureQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_h2oQuantity) ! For hydrostatic
            h2oVectorIndex=decoration(decoration(subtree(1,gson)))
            h2oQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_geocAltitudeQuantity) ! For hydrostatic
            geocAltitudeVectorIndex=decoration(decoration(subtree(1,gson)))
            geocAltitudeQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_refGPHQuantity) ! For hydrostatic
            refGPHVectorIndex=decoration(decoration(subtree(1,gson)))
            refGPHQuantityIndex=decoration(decoration(decoration(subtree(2,gson))))
          case (f_explicitValues) ! For explicit fill
            valuesNode=subtree(j,key)
          case (f_sourceGrid)
            gridIndex=decoration(decoration(gson))
          case (f_maxIterations)      ! For hydrostatic fill
            call expr(subtree(2,subtree(j,key)), unitAsArray,valueAsArray)
            if (all(unitAsArray(1) /= (/PHYQ_Dimensionless,PHYQ_Invalid/))) &
              & call Announce_error( key, badUnitsForMaxIterations)
            maxIterations=valueAsArray(1)
          case (f_spread) ! For explicit fill, note that gson here is not same as others
            if (node_id(gson) == n_set_one) then
              spread=.TRUE.
            else
              spread=decoration(subtree(2,gson)) == l_true
            endif
          end select
        end do                  ! Loop over arguments to fill instruction

        ! Now call various routines to do the filling
        quantity=>GetVectorQtyByTemplateIndex(vectors(vectorIndex),quantityIndex)
        select case (fillMethod)


        case (l_hydrostatic) ! -------------------------- Hydrostatic fills ------
          ! Need a temperature and a refgph quantity
          if (.not.all(got( (/ f_refGPHQuantity, f_temperatureQuantity /)))) &
            call Announce_Error(key,needTempREFGPH)

          temperatureQuantity => GetVectorQtyByTemplateIndex( &
            &  vectors(temperatureVectorIndex), temperatureQuantityIndex)
          if (temperatureQuantity%template%quantityType /= l_Temperature) &
            & call Announce_Error (key, badTemperatureQuantity)

          refGPHQuantity => GetVectorQtyByTemplateIndex( &
            & vectors(refGPHVectorIndex), refGPHQuantityIndex)
          if (refGPHQuantity%template%quantityType /= l_refGPH) &
            & call Announce_Error (key, badrefGPHQuantity)

          if (quantity%template%quantityType==l_ptan) then
            if (.not. got(f_geocAltitudeQuantity)) &
              & call Announce_Error( key, needGeocAltitude )
            geocAltitudeQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(geocAltitudeVectorIndex), geocAltitudeQuantityIndex)
            if (geocAltitudeQuantity%template%quantityType /= l_tngtgeocAlt) &
              & call Announce_Error( key, badGeocAltitudeQuantity )
            if (.not. got(f_h2oQuantity)) &
              & call Announce_Error( key, needH2O )
            h2oQuantity => GetVectorQtyByTemplateIndex( &
              & vectors(h2oVectorIndex), h2oQuantityIndex)
            if (.not. ValidateVectorQuantity(h2oQuantity, &
              & quantityType=(/l_vmr/), molecule=(/l_h2o/)))&
              & call Announce_Error( key, badGeocAltitudeQuantity )
          else
            geocAltitudeQuantity=>NULL()
            h2oQuantity=>NULL()
          endif
          call FillVectorQtyHydrostatically(key, quantity, temperatureQuantity, &
            & refGPHQuantity, h2oQuantity, geocAltitudeQuantity, maxIterations)          

        case (l_special) ! ------------------ Special fills for some quantities --
          select case (quantity%template%quantityType)
          case (l_losVel)
            if (.not. any(got( (/f_tngtECI, f_scECI, f_scVel/) ))) then
              call Announce_error(key, badlosVelFill)
            else
              tngtECIQuantity=> GetVectorQtyByTemplateIndex( &
                & vectors(tngtECIVectorIndex), tngtECIQuantityIndex)
              scECIQuantity=> GetVectorQtyByTemplateIndex( &
                & vectors(scECIVectorIndex), scECIQuantityIndex)
              scVelQuantity=> GetVectorQtyByTemplateIndex( &
                & vectors(scVelVectorIndex), scVelQuantityIndex)
              call FillLOSVelocity(key, quantity, tngtECIQuantity, &
                & scECIquantity, scVelQuantity)
            endif
          case default
            call Announce_error(key, noSpecialFill)
          end select

        case (l_gridded) ! --------------------- Fill from gridded data --
          if (.not. got(f_sourceGrid)) call Announce_Error(key,noSourceGridGiven)
          call FillVectorQuantityFromGrid(quantity,griddedData(gridIndex),errorCode)
          if (errorCode/=0) call Announce_error(key,errorCode)

        case (l_l2gp) ! ---------------------- Fill from L2GP quantity ---
          if (.NOT. got(f_sourceL2GP)) call Announce_Error(key,noSourceL2GPGiven)
          call FillVectorQuantityFromL2GP(quantity,l2gpDatabase(l2gpIndex),errorCode)
          if (errorCode/=0) call Announce_error(key,errorCode)


        case (l_l2aux) ! --------------------- Fill from L2AUX quantity --
          if (.NOT. got(f_sourceL2AUX)) call Announce_Error(key,noSourceL2AUXGiven)
!          call FillVectorQuantityFromL2AUX(quantity,l2auxDatabase(l2auxIndex),errorCode)
          if (errorCode/=0) call Announce_error(key,errorCode)


        case (l_explicit) ! -------------------- Explicity fill from l2cf -
          if (.not. got(f_explicitValues)) call Announce_Error(key, &
            & noExplicitValuesGiven)
          call ExplicitFillVectorQuantity ( quantity, valuesNode, spread )


        case (l_l1b)                    ! Fill from L1B data
          call FillVectorQuantityFromL1B ( key, quantity, chunks(chunkNo), l1bInfo )


        case default
!          call MLSMessage(MLSMSG_Error,ModuleName,'This fill method not yet implemented')
          call Announce_error(key,0, &
			 & 'This fill method not yet implemented')
        end select
        
        ! End of fill operations

      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if

        case ( s_snoop )
          call Snoop( key=key, vectorDatabase=vectors)
      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

    if (ERROR/=0) then
	 !	call MLSMessage(MLSMSG_Error,ModuleName,'Problem with Fill section')
          call Announce_error(key,0, &
			 & 'Problem with Fill section')
	endif

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 ) then
        call dump ( vectors )
      end if
      call trace_end ( "MLSL2Fill" )
    end if
    if ( timing ) call sayTime

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for MLSL2Fill =" )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Fill

  ! =====     Private Procedures     =====================================

  !=============================== FillVector ==========================
  subroutine FillVector(errorCode, inArray, Vector, arrayType, vectorType, NumQtys, &
    & numChans, numSurfs, numInstances, qtiesStart, dim_order)
    !=============================== FillVector ==========================

    ! Fill the vector Vector with values taken from the array inArray
    ! in a manner that depends on their respective types:
!
    !        arrayType        vectorType            operation
    !          l2gp              l2gp       vector(:,:,:) = inArray(:,:,:)

    ! It is assumed that the rank3 inArray is filled from
    ! inArray(1, 1, 1) to inArray(numChans, numSurfs, numInstances)
!
    ! With the redefinition of Vector%quantities to be a rank2 object
    ! we are faced with a problem: how to fill a rank 2 object from a rank 3 one?
    ! If the rank 3 object is full, then it use the following trick:
    ! Vector(1:numChans*numSurfs, 1:numInstances) = 
    !                   inArray(1:numChans, 1:numSurfs, 1:numInstances)

    character (Len=*), intent(IN) ::        arrayType, vectorType
    integer, intent(OUT) ::                 errorCode		! zero unless an error
    ! REAL(r8), POINTER, DIMENSION(:,:,:) ::  inPointer
    real(r8), dimension(:,:,:) ::           inArray
    type(Vector_T), intent(OUT) ::          Vector
    integer, intent(IN) ::                  NumQtys, qtiesStart
    integer, intent(IN) ::                  numInstances, numSurfs, numChans
    integer, intent(IN), optional, dimension(:) :: dim_order

    ! Private
    character (LEN=LineLen) ::              msr
    integer ::                              qty
    character (LEN=4) ::                    qtyChar, errorCodeChar

    ! Sanity checks:
    if(numChans.eq.0) then
      !	call announce_error ( ErrorInFillVector, numChansisZero )
      errorCode = numChansisZero
      return
    elseif(numSurfs.eq.0) then
      !	call announce_error ( ErrorInFillVector, numSurfsisZero )
      errorCode = numSurfsisZero
      return
    elseif(numInstances.eq.0) then
      !	call announce_error ( ErrorInFillVector, numInstancesisZero )
      errorCode = numInstancesisZero
      return
    endif

    Vector%template%NoQuantities = NumQtys
    select case (arrayType(:4) // VectorType(:4))
    case ('l2gpl2gp')
      do qty = qtiesStart, qtiesStart - 1 + NumQtys
        write(qtyChar, '(I4)') qty
        Vector%quantities(qty)%template%noChans = numChans
        Vector%quantities(qty)%template%noSurfs = numSurfs
        Vector%quantities(qty)%template%noInstances = numInstances    
        !     CALL put_3d_view(inPointer, &
        !        & numChans, &
        !        & numSurfs, &
        !        & numInstances, &
        !        & Vector%quantities(qty)%values)
        !     Vector%quantities(qty)%values = inPointer
        !	IF(numChans.EQ.1) THEN
        !        	Vector%quantities(qty)%values = inPointer(:, :, 1)
        !	ELSEIF(numSurfs.EQ.1) THEN
        !        	Vector%quantities(qty)%values = inPointer(:, 1, :)
        !	ELSEIF(numInstances.EQ.1) THEN
        !        	Vector%quantities(qty)%values = inPointer(1, :, :)
        !        ELSE
        !		call announce_error ( ErrorInFillVector, objIsFullRank3 )
        !       		RETURN
        !	ENDIF
	call squeeze(errorCode, inArray, Vector%quantities(qty)%values)
        !     WRITE(errorCodeChar, '(I4)') errorCode
        !     msr = 'Error #' // errorCodeChar // ' in squeezing l2gp into quantity ' // qtyChar
        !     IF(errorCode /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
        !        & msr)
      enddo
    case ('l2aul2au')
      do qty = qtiesStart, qtiesStart - 1 + NumQtys
        write(qtyChar, '(I4)') qty
        Vector%quantities(qty)%template%noChans = numChans
        Vector%quantities(qty)%template%noSurfs = numSurfs
        Vector%quantities(qty)%template%noInstances = numInstances
        if(present(dim_order)) then
          call squeeze(errorCode, inArray, Vector%quantities(qty)%values, dim_order)
        else    
          call squeeze(errorCode, inArray, Vector%quantities(qty)%values)
        endif
        !     WRITE(errorCodeChar, '(I4)') errorCode
        !     msr = 'Error #' // errorCodeChar // ' in squeezing l2aux into quantity ' // qtyChar
        !     IF(errorCode /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
        !        & msr)
      enddo
    case default
      ! FillVector not yet written to handle these cases
    end select

  end subroutine FillVector

  !=============================== FillOL2GPVector ==========================
  subroutine FillVectorQuantityFromGrid(quantity,grid, errorCode)
    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (GriddedData_T), intent(in) :: GRID ! Grid to fill it from
    integer, intent(out) :: ERRORCODE   ! Error code (one of constants defined above)

    ! Local variables

  end subroutine FillVectorQuantityFromGrid

  !=============================== FillOL2GPVector ==========================
  subroutine FillVectorQuantityFromL2GP(quantity,l2gp, errorCode)

    ! If the times, pressures, and geolocations match, fill the quantity with
    ! the appropriate subset of profiles from the l2gp

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (L2GPData_T), intent(in) :: L2GP ! L2GP to fill from
    integer, intent(out) :: errorCode ! Error code

    ! Local parameters
    real(r8), parameter:: TOLERANCE=0.05 ! Tolerence for angles

    ! Local variables
    integer ::    FIRSTPROFILE, LASTPROFILE
    integer, dimension(1) :: FIRSTPROFILEASARRAY

    errorCode=0
    ! Make sure this quantity is appropriate
    if (.not. ValidateVectorQuantity(quantity, coherent=.TRUE., stacked=.TRUE., &
      & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if ( (quantity%template%noChans/=l2gp%nFreqs) .and. &
      &  ((quantity%template%noChans/=1) .or. (l2gp%nFreqs/=0)) ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if ( quantity%template%noSurfs /= l2gp%nLevels ) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if ( quantity%template%verticalCoordinate == l_pressure) then
      if ( any(ABS(-LOG10(quantity%template%surfs(:,1))+ &
        & LOG10(l2gp%pressures)) > TOLERANCE)) then
        errorCode=vectorWontMatchL2GP
        return
      end if
    else                                ! Must be l_zeta
      if ( any(ABS(quantity%template%surfs(:,1)+ &
        & LOG10(l2gp%pressures)) > TOLERANCE)) then
        errorCode=vectorWontMatchL2GP
        return
      end if
    end if

    ! Attempt to match up the first location
    firstProfileAsArray=MINLOC(ABS(quantity%template%phi(1,1)-l2gp%geodAngle))
    firstProfile=firstProfileAsArray(1)

    ! Well, the last profile has to be noInstances later, check this would be OK
    lastProfile=firstProfile+quantity%template%noInstances-1
    if (lastProfile > l2gp%nTimes) then
      errorCode=vectorWontMatchL2GP
      return
    endif

    ! Now check that geodAngle's are a sufficient match
    if (any(abs(l2gp%geodAngle(firstProfile:lastProfile)-&
      &         quantity%template%phi(1,:)) > tolerance)) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    if (any(abs(l2gp%time(firstProfile:lastProfile)- &
      &         quantity%template%time(1,:)) > tolerance)) then
      errorCode=vectorWontMatchL2GP
      return
    end if

    quantity%values=RESHAPE(l2gp%l2gpValue(:,:,firstProfile:lastProfile),&
      & (/quantity%template%noChans*quantity%template%noSurfs,&
      &   quantity%template%noInstances/))

  end subroutine FillVectorQuantityFromL2GP

  ! ------------------------------------------- FillLOSVelocity ---
  subroutine FillLOSVelocity ( key, qty, tngtECI, scECI, scVel)
    ! A special fill from geometry arguments
    integer, intent(in) :: KEY
    type (VectorValue_T), intent(inout) :: QTY
    type (VectorValue_T), intent(in) :: TNGTECI
    type (VectorValue_T), intent(in) :: SCECI
    type (VectorValue_T), intent(in) :: SCVEL


    ! Local variables
    integer :: MAF                      ! Loop counter
    integer :: MIF                      ! Loop counter
    integer :: noMAFs                   ! Number of major frames
    integer :: noMIFs                   ! Number of minor frames for this module
    integer :: x,y,z                    ! Indicies into the vectors

    real (r8), dimension(3) :: tngtVel   ! Due to rotation of earth
    real (r8), dimension(3) :: los       ! Normalised line of sight vector

    ! Executable code
    ! First check that things are OK.
    if ( (qty%template%quantityType /= l_losVel) .or. &
      &  (tngtECI%template%quantityType /= l_tngtECI) .or. &
      &  (scECI%template%quantityType /= l_scECI) .or. &
      &  (scVel%template%quantityType /= l_scVel)) then
      call Announce_Error(key, badLOSVelFill)
      return
    endif

    if ( qty%template%instrumentModule /= tngtECI%template%instrumentModule ) then
      call Announce_Error(key, badLOSVelFill)
      return
    endif

    noMAFs = qty%template%noInstances
    noMIFs = qty%template%noSurfs

    do maf = 1, noMAFs
      do mif = 1, noMIFs

        ! First compute the tangent point velocity in ECI coordinates due 
        ! to the rotation of the earth.  This no doubt makes approximations
        ! due to the slight non alignment between the earth's rotation axis and
        ! the ECI z axis, but I'm going to ignore this.

        ! Work out the indices in 3*mif,maf space
        x = 1 + 3*(mif-1)
        y = x+1
        z = x+2

        tngtVel= omega* (/ -tngtECI%values(y,maf), &
          &                 tngtECI%values(x,maf), 0.0_r8 /)

        ! Now compute the line of sight direction normal
        los = tngtECI%values(x:z,maf) - scECI%values(x:z,maf)
        los = los / sqrt(sum(los**2))

        ! Now compute the net velocity in this direction.  For the moment I'll
        ! assume +ve means the sc and tp are moving apart, and -ve that they're
        ! getting closer.

        qty%values(mif,maf) = dot_product(tngtVel, los) - &
          &                   dot_product(scVel%values(x:z,maf), los)

        ! Note that even though x,y,z have been used up to now for a GHz/THz
        ! minor frame quantity, they're OK with this sc one too.
      end do
    end do
  end subroutine FillLOSVelocity

  ! ------------------------------------- FillVectorHydrostatically ----
  subroutine FillVectorQtyHydrostatically(key, quantity, &
    & temperatureQuantity, refGPHQuantity, h2oQuantity, &
    & geocAltitudeQuantity, maxIterations)
    ! Various hydrostatic fill operations
    integer, intent(in) :: key          ! For messages
    type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
    type (VectorValue_T), intent(in) :: TEMPERATUREQUANTITY
    type (VectorValue_T), intent(in) :: REFGPHQUANTITY
    type (VectorValue_T), pointer :: H2OQUANTITY
    type (VectorValue_T), pointer :: GEOCALTITUDEQUANTITY
    integer, intent(in) :: MAXITERATIONS
    ! H2OQuantity and GeocAltitudeQuantity have to be pointers
    ! as they may be absent.

    ! Local variables

    ! Executable code

    if ( toggle(gen) ) call trace_begin ( "FillVectorQtyHydrostatically", key )

    select case (quantity%template%quantityType)
    case (l_gph)
      if ( (temperatureQuantity%template%noSurfs /= &
        &   quantity%template%noSurfs) .or. &
        &  (refGPHQuantity%template%noInstances /= &
        &   quantity%template%noInstances) .or. &
        &  (temperatureQuantity%template%noInstances /= &
        &   quantity%template%noInstances) ) then
        call Announce_Error ( key, nonConformingHydrostatic )
		    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      endif
      if ((any(quantity%template%surfs /= temperatureQuantity%template%surfs)) .or. &
        & (any(quantity%template%phi /= temperatureQuantity%template%phi)) .or. &
        & (any(quantity%template%phi /= refGPHQuantity%template%phi)) ) then
        call Announce_Error ( key, nonConformingHydrostatic )
		    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      endif
      call GetBasisGPH(temperatureQuantity, refGPHQuantity, quantity%values)
    case (l_ptan)
      if ( (quantity%template%noInstances /= &
        &   temperatureQuantity%template%noInstances) .or. &
        &  (quantity%template%noInstances /= &
        &   refGPHquantity%template%noInstances) .or. &
        &  (quantity%template%noInstances /= &
        &   h2oQuantity%template%noInstances) ) then
        call Announce_Error ( key, nonConformingHydrostatic )
		    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      endif
      if ((any(refGPHquantity%template%phi /= temperatureQuantity%template%phi)) .or. &
        & (any(h2oQuantity%template%phi /= temperatureQuantity%template%phi)) ) then
        call Announce_Error ( key, nonConformingHydrostatic )
		    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      endif
      if ( (.not. ValidateVectorQuantity(quantity, minorFrame=.true.) ) .or. &
        &  (.not. ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ) .or. &
        &  (quantity%template%instrumentModule /= &
        &   geocAltitudeQuantity%template%instrumentModule) )  then
        call Announce_Error (key, nonConformingHydrostatic )
		    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically")
        return
      end if
      call GetHydrostaticTangentPressure(quantity, temperatureQuantity,&
        & refGPHQuantity, h2oQuantity, geocAltitudeQuantity, maxIterations)
    case default
!      call MLSMessage(MLSMSG_Error, ModuleName, 'No such fill yet')
          call Announce_error(0, 0, &
			 & 'No such fill yet')
    end select

    if ( toggle(gen) ) call trace_end ( "FillVectorQtyHydrostatically" )

  end subroutine FillVectorQtyHydrostatically

  !=============================== FillPrevDefd ==========================
  subroutine FillPrevDefd(OldVector, PrevDefdQt, qtyTemplates, &
    & Output, QuantityName, qtiesStart, chunkNo, errorCode)
    !=============================== FillPrevDefd ==========================

    ! If the times, pressures, and geolocations match,
    ! fill the vector Output with values taken from the appropriate quantity in
    ! previously defined vector OldVector
    ! 

    integer, intent(IN) ::                            PrevDefdQt
    integer, intent(IN) ::                            qtiesStart
    type(Vector_T), intent(IN) ::                     OldVector
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
    type(Vector_T), intent(INOUT) ::                  Output
    character*(*), intent(IN) ::                      QuantityName
    ! If done chunk-by-chunk, the following is the chunk number
    ! Otherwise, chunkNo should = -1
    integer, intent(in) ::                        chunkNo
    integer, intent(OUT) ::                            errorCode	! if error

    ! Local variables
    !::::::::::::::::::::::::: LOCALS :::::::::::::::::::::

    real(r8) ::                                       TOLERANCE
    parameter(TOLERANCE=0.05)

    ! A wrong version compared ChunkNo to OldVector%ChunkNumbers
    logical ChunkNumberIsTrustworthy
    parameter(ChunkNumberIsTrustworthy=.false.)

    type (QuantityTemplate_T) ::                      PDQTemplate	! Prev. def'd Quantity Template
    type (QuantityTemplate_T) ::                      OQTemplate	! Output Quantity Template
    integer ::                                        i
    integer ::                                        firstProfile, lastProfile
    integer ::                                        noL2GPValues=1
    logical ::                                        TheyMatch
    real(r8) ::                                       phiMin, phiMax
    real(r8) ::                                       phi_TOLERANCE

    PDQTemplate = qtyTemplates(OldVector%TEMPLATE%QUANTITIES(PrevDefdQt))
    OQTemplate = qtyTemplates(Output%TEMPLATE%QUANTITIES(qtiesStart))

    if(OQTemplate%noInstances <= 0) then
      errorCode=zeroProfilesFound
      return
    elseif(PDQTemplate%noInstances <= 0) then
      errorCode=zeroProfilesFound
      return
    endif
    ! Chunk-by-chunk, or all chunks at once?
    if(chunkNo == -1) then
      firstProfile = 1
      lastProfile = PDQTemplate%noInstances
    elseif(ChunkNumberIsTrustworthy) then
      !	lastProfile = 1
      !   firstProfile = PDQTemplate%noInstances
      !	DO i=1, PDQTemplate%noInstances
      !        	IF(chunkNo == OldVector%chunkNumber(i)) THEN
      !				firstProfile = MIN(firstProfile, i)
      !				lastProfile = MAX(lastProfile, i)
      !         ENDIF
      !   ENDDO
      !Can't get here unless ChunkNumberIsTrustworthy has been reset
      call announce_error(0, miscellaneous_err, &
	& "Programming error in Module Fill, SUBROUTINE FillPrevDefd " &
	& // "ChunkNumberIsTrustworthy must be FALSE")
    else
      ! Instead of comparing OldVector%chunkNumbers to ChunkNo
      ! we will compare geodetic angles to phi
      phiMin = OQTemplate%phi(1, 1)
      phiMax = OQTemplate%phi(1, 1)
      do i=1, OQTemplate%noInstances
        phiMin = min(phiMin, OQTemplate%phi(1, i))
        phiMax = max(phiMax, OQTemplate%phi(1, i))
      enddo
      phi_TOLERANCE = TOLERANCE *(phiMax-phiMin) / OQTemplate%noInstances
      if(phi_TOLERANCE <= 0.D0) then
        errorCode=zeroGeodSpan
        return
      endif
      lastProfile = 1
      firstProfile = PDQTemplate%noInstances
      do i=1, PDQTemplate%noInstances
        if(phiMin <= (PDQTemplate%phi(1, i) + phi_TOLERANCE)) then
          firstProfile = min(firstProfile, i)
        endif
        if(phiMax >= (PDQTemplate%phi(1, i) - phi_TOLERANCE)) then
          lastProfile = max(lastProfile, i)
        endif
      enddo
    endif
!
    TheyMatch = OQTemplate%noInstances .eq. (lastProfile-firstProfile+1)
    TheyMatch = TheyMatch &
      & .and. &
      & (OQTemplate%stacked .eqv. PDQTemplate%stacked) &
      & .and. &
      & (OQTemplate%coherent .eqv. PDQTemplate%coherent) &
      & .and. &
      & (OQTemplate%regular .eqv. PDQTemplate%regular) 
    if(TheyMatch) then
      do i = firstProfile, lastProfile
        if( &
          &   nearBy(PDQTemplate%geodLat(1,i), OQTemplate%geodLat(1,i)) &
          & .and. &
          &   nearBy(PDQTemplate%lon(1,i), OQTemplate%lon(1,i)) &
          & .and. &
          &   nearBy(PDQTemplate%solarTime(1,i), OQTemplate%solarTime(1,i)) &
          &) then
        else
          TheyMatch = .false.
          exit
        endif
      end do
    endif
    if( TheyMatch &
      & .and. &
      & PDQTemplate%noSurfs == OQTemplate%noSurfs &
      & .and. &
      & PDQTemplate%noSurfs == 1 ) then
      do i = 1, OQTemplate%noSurfs
        if( &
          &   nearBy(PDQTemplate%surfs(i,1), OQTemplate%surfs(i,1)) &
          & ) then
        else
          TheyMatch = .false.
          exit
        endif
      end do
    endif
    if(TheyMatch) then
      Output%quantities(qtiesStart)%values = OldVector%quantities(PrevDefdQt)%values
    else
      errorCode=vectorWontMatchPDef
    endif
  end subroutine FillPrevDefd

  !=============================== FillOL2AUXVector ==========================
  subroutine FillOL2AUXVector(OldL2AUXData, qtyTemplates, &
    & Output, QuantityName, qtiesStart, chunkNo, errorCode)
    !=============================== FillOL2AUXVector ==========================

    ! If the times, pressures, and geolocations match,
    ! fill the vector Output with values taken from the
    ! Old L2AUX vector OldL2AUXData
    ! 

    integer, intent(IN) ::                            qtiesStart
    type(L2AUXData_T) ::                               OldL2AUXData
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
    !TYPE(L2GPData_T), DIMENSION(:), POINTER ::        L2GPDatabase
    type(Vector_T), intent(INOUT) ::                  Output
    character*(*), intent(IN) ::                      QuantityName
    integer, intent(in) :: chunkNo
    integer, intent(OUT) ::                            errorCode	! if error

    ! Local variables
    !::::::::::::::::::::::::: LOCALS :::::::::::::::::::::
    type (QuantityTemplate_T) ::                      OQTemplate	! Output Quantity Template
    integer ::                                        OQType	! Ouptut quantity type
    integer ::                                        OQNVals	! No. Ouptut quantity vals
    integer ::                                        ChanDim	! which dim no. is channel
    integer ::                                        IntFreqDim	! which dim no. is int. frq.
    integer ::                                        USBDim	! which dim no. is USB
    integer ::                                        LSBDim	! which dim no. is LSB
    integer ::                                        MAFDim	! which dim no. is MAF
    integer ::                                        MIFDim	! which dim no. is MIF
    integer ::                                        L2NVals	! No. l2 aux vals
    integer ::                                        i
    logical ::                                        TheyMatch
    ! This will let us re-order the vector dimensions differently from the l2aux
    integer, dimension(3) ::                          dim_order

    ! Properties of Output Vector
    OQTemplate = qtyTemplates(Output%TEMPLATE%QUANTITIES(qtiesStart))
    OQType = OQtemplate%quantityType
    OQNVals = OQtemplate%noInstances*OQtemplate%noSurfs*OQtemplate%noChans

    ! Properties of l2 aux
    ChanDim = 0
    IntFreqDim = 0
    USBDim = 0
    LSBDim = 0
    MAFDim = 0
    MIFDim = 0
    L2NVals = 1
    do i=1, L2AUXRank
      L2NVals = L2NVals*OldL2AUXData%dimensions(i)%noValues
      select case(OldL2AUXData%dimensions(i)%dimensionFamily)
      case(L_Channel)
        ChanDim = i
        dim_order(1) = i
      case(L_IntermediateFrequency)
        IntFreqDim = i
      case(L_USBFrequency)
        USBDim = i
      case(L_LSBFrequency)
        LSBDim = i
      case(L_MIF)
        MIFDim = i
      case(L_MAF)
        MAFDim = i
      case DEFAULT	! We are not yet interested in these dimensions
      end select

    enddo

    ! Check that the dimensions match up and are of the same family
    TheyMatch = OQNVals == L2NVals
    if(OQTemplate%minorFrame) then
      TheyMatch = TheyMatch .and. MIFDim > 0
    endif
    if(OQTemplate%noChans > 0) then
      TheyMatch = TheyMatch .and. &
        & (  ChanDim > 0   )
    endif

    if(TheyMatch) then
      call FillVector(errorCode, OldL2AUXData%values, Output, &
        & 'l2aux', 'l2aux', 1, &
        & max(OldL2AUXData%dimensions(1)%noValues, 1), &
        & max(OldL2AUXData%dimensions(2)%noValues, 1), &
        & max(OldL2AUXData%dimensions(3)%noValues, 1), &
        & qtiesStart)
    else
      !   CALL MLSMessage(MLSMSG_Error, ModuleName, &
      !        & 'Vector and old L2AUX do not match')
      errorCode=cantFillFromL2AUX
    endif

  end subroutine FillOL2AUXVector

  !=============================================== ExplicitFillVectorQuantity ==
  subroutine ExplicitFillVectorQuantity(quantity, valuesNode, spread)

    ! This routine is called from MLSL2Fill to fill values from an explicit
    ! fill command line

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
    integer, intent(in) :: VALUESNODE ! Tree node
    logical, intent(in) :: SPREAD ! One instance given, spread to all

    ! Local variables
    integer :: K                ! Loop counter
    integer, DIMENSION(2) :: unitAsArray ! Unit for value given
    real (r8), DIMENSION(2) :: valueAsArray ! Value given
    
    ! Executable code

    if (spread) then      ! 1 instance given, spread to all instances

      ! Check we have the right number of values
      if ( (nsons(valuesNode)-1 /= quantity%template%instanceLen) .or. &
        &  (.not. quantity%template%regular)) &
        & call Announce_error(valuesNode,invalidExplicitFill)

      ! Loop over the values
      do k=1,nsons(valuesNode)-1
        ! Get value from tree
        call expr(subtree(k+1,valuesNode),unitAsArray,valueAsArray)
        ! Check unit OK
       if ( (unitAsArray(1) /= quantity%template%unit) .and. &
          &  (unitAsArray(1) /= PHYQ_Dimensionless) ) &
          & call Announce_error(valuesNode,badUnitsForExplicit)
        ! Store value
        quantity%values(k,:)=valueAsArray(1)
      end do

    else                  ! Not spread, fill all values

      ! Check we have the right number of values
      if (nsons(valuesNode)-1 /= &
        & quantity%template%noInstances*quantity%template%instanceLen) &
        & call Announce_error(valuesNode,invalidExplicitFill)

      ! Loop over values
      do k=1,nsons(valuesNode)-1
        ! Get value from tree
        call expr(subtree(k+1,valuesNode),unitAsArray,valueAsArray)
        ! Check unit OK
        if ( (unitAsArray(1) /= quantity%template%unit) .and. &
          &  (unitAsArray(1) /= PHYQ_Dimensionless) ) &
          & call Announce_error(valuesNode,badUnitsForExplicit)
        ! Store value
        quantity%values(mod(k-1,quantity%template%instanceLen)+1,&
          &             (k-1)/quantity%template%instanceLen+1)=&
          & valueAsArray(1)
      end do
    endif
  end subroutine ExplicitFillVectorQuantity

  ! ----------------------------------------- FillVectorQuantityFromL1B ----
  subroutine FillVectorQuantityFromL1B ( root, quantity, chunk, l1bInfo )
    integer, intent(in) :: root
    type (VectorValue_T), INTENT(INOUT) :: QUANTITY
    type (MLSChunk_T), INTENT(IN) :: CHUNK
    type (l1bInfo_T), INTENT(IN) :: L1BINFO

    ! Local variables
    character (len=80) :: NAMESTRING
    integer :: fileID, FLAG, NOMAFS
    type (l1bData_T) :: L1BDATA

    ! Executable code

    if ( toggle(gen) ) call trace_begin ("FillVectorQuantityFromL1B",root)

    fileID=l1bInfo%l1bOAID
    select case ( quantity%template%quantityType )
    case ( l_radiance )
      call GetSignalName ( quantity%template%signal, nameString, noChannels=.TRUE. )
      fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString )
    case ( l_tngtECI )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpECI'
    case ( l_tngtGeodAlt )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpGeodAlt'
    case ( l_tngtGeocAlt )
      call GetModuleName( quantity%template%instrumentModule,nameString )
      nameString=TRIM(nameString)//'.tpGeocAlt'
    case ( l_scECI )
      nameString='scECI'
    case ( l_scVel )
      nameString='scVel'
    case ( l_scGeocAlt )
      nameString='scGeocAlt'
    case default
      call Announce_Error(root, cantFillFromL1B)
    end select

    call ReadL1BData ( fileID , nameString, l1bData, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    ! We'll have to think about `bad' values here .....
    if ( flag /= 0 ) then
      call Announce_Error(errorReadingL1B, root)
		    if ( toggle(gen) ) call trace_end ( "FillVectorQuantityFromL1B")
      return
    end if
    quantity%values = RESHAPE(l1bData%dpField, &
      & (/ quantity%template%instanceLen, quantity%template%noInstances /) )
    call DeallocateL1BData(l1bData,flag)

    if (toggle(gen) ) call trace_end( "FillVectorQuantityFromL1B" )
  end subroutine FillVectorQuantityFromL1B

  !=============================== nearby ==========================
  function nearby(x, y, inRelSmall, inOverallTol)
    !=============================== nearby ==========================
    ! This functions returns TRUE if x and y are "nearby" comparing either
    ! their relative difference (proportionally small) 
    ! or their differences with an overall tolerance
    real(r8), intent(IN) ::              x, y
    real(r8), intent(IN), optional ::    inRelSmall
    real(r8), intent(IN), optional ::    inOverallTol

    ! result
    logical ::                           nearby

    ! Private
    real(r8) ::                          small=1.D-32
    real(r8) ::                          tolerance=0.D0

    if(present(inRelSmall)) then
      small = inRelSmall
    endif

    if(present(inOverallTol)) then
      tolerance = inOverallTol
    endif

    if(abs(x-y) <= tolerance) then
      nearby = .true.
    elseif(x*y /= 0) then
      nearby = abs(x-y) .lt. small*sqrt(abs(x))*sqrt(abs(y))
    elseif(x /= 0) then
      nearby = abs(x-y) .lt. small*abs(x)
    else
      nearby=.true.
    endif

  end function nearby

  !=============================== whatQuantityNumber ==========================
  function whatQuantityNumber(x, y, xVectorTemplate)
    !=============================== whatQuantityNumber ==========================
    ! This functions returns the quantity number of quantity y in vector x
    ! where x and y are sub-rosa indexes of their actual names;
    ! i.e., they are specified in the cf as name[x].name[y]
    ! where name[x] = get_char(x), name[y] = get_char[y]
!
    ! If the quantity is not found, it returns FAILED (-1)
    integer, intent(IN) ::               x, y
    type (VectorTemplate_T)  ::          xVectorTemplate

    ! result
    integer ::                           whatQuantityNumber

    ! Private
    integer, parameter ::                FAILED=-1
    integer              ::              qty

    whatQuantityNumber = FAILED
    do qty = 1, xVectorTemplate%NoQuantities
      if(xVectorTemplate%Name == y) then
        whatQuantityNumber = qty
        exit
      endif
    enddo

  end function whatQuantityNumber

  !=============================== squeeze ==========================
  subroutine squeeze(errorCode, source, sink, source_order)
    !=============================== squeeze ==========================
    ! takes a rank 3 object source and returns a rank2 object sink
    ! source(1..n1, 1..n2, 1..n3) -> sink(1..n1*n2, 1..n3)
    ! unless it can't--then it returns errorCode /= 0
    ! One reason it may fail: shape of sink too small
!
    ! Assuming that shape(source) = {n1, n2, n3}
    !     =>          shape(sink) = {m1, m2}
    ! then we must further assume (else set errorCode)
    ! n1*n2 <= m1
    ! n3 <= m2
!
    ! if source_order is present, it is a permutation of (1 2 3)
    ! such that before squeezing, the source is re-ordered:
    ! temp(1, 2, 3) = source(order(1), order(2), order(3))
    ! and then
    ! temp(1..n1, 1..n2, 1..n3) -> sink(1..n1*n2, 1..n3)
!
    ! (A future improvement might take as optional arguments
    !  integer arrays source_shape, sink_shape, 
    !  or else shape-params n1, n2, m1)
    !--------Argument--------!
    real(r8), dimension(:,:,:), intent(IN) :: source
    real(r8), dimension(:,:), intent(OUT)  :: sink
    integer, intent(OUT)                   :: errorCode
    integer, intent(IN), optional, dimension(:) :: source_order

    !----------Local vars----------!

    integer, dimension(4) :: source_shape
    integer, dimension(4) :: sink_shape
    integer, dimension(4) :: temp_shape
    integer::i,icode,offset
    real(r8), dimension(:,:,:), allocatable :: temp
    !----------Executable part----------!
    source_shape(1:3) = shape(source)
    sink_shape(1:2) = shape(sink)

    if (source_shape(1) == 0) then
      errorCode = n1_is_zero
      return
    elseif (source_shape(2) == 0) then
      errorCode = n2_is_zero
      return
    elseif (source_shape(3) == 0) then
      errorCode = n3_is_zero
      return
    elseif (sink_shape(1) < source_shape(1)*source_shape(2)) then
      errorCode = m1_too_small
      return
    elseif (sink_shape(2) < source_shape(3)) then
      errorCode = m2_too_small
      return
    else
      errorCode = 0
    endif

    if(present(source_order)) then
      !        Check that source_order is a legal permutation of (1 2 3)
      !        using trick: sum and product must each equal 6
      if(source_order(1)+source_order(2)+source_order(3) /= 6 &
        & .or. &
        & source_order(1)*source_order(2)*source_order(3) /= 6 ) then
        errorCode=not_permutation
        return
      endif
      do i=1, 3
        temp_shape(i) = source_shape(source_order(i))
      enddo
      allocate(temp(temp_shape(1), temp_shape(2), temp_shape(3)), &
        & STAT=errorCode)
      if(errorCode /= 0) then
        errorCode=allocation_err
        return
      endif
      temp = reshape(source, temp_shape(1:3), order=source_order(1:3))
      sink = reshape(temp, sink_shape(1:2))
      deallocate(temp, &
        & STAT=errorCode)
      if(errorCode /= 0) then
        errorCode=deallocation_err
        return
      endif
    else
      sink = reshape(source, sink_shape(1:2))
    endif

  end subroutine squeeze

!
!
  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( where, CODE , ExtraMessage)
    integer, intent(in) :: where   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    character (LEN=*), optional :: ExtraMessage

    error = max(error,1)
    call output ( '***** At ' )

	if(where > 0) then
    call print_source ( source_ref(where) )
		else
    call output ( '(no lcf tree available)' )
		endif

    call output ( ': ' )
    call output ( "The " );

	if(where > 0) then
    call dump_tree_node ( where, 0 )
		else
    call output ( '(no lcf node available)' )
		endif

    select case ( code )
    case ( badUnitsForExplicit )
      call output ( " has inappropriate units for Fill instruction.", advance='yes' )
    case ( wrong_number )
      call output ( " command does not have exactly one field.", advance='yes' )
    case ( unknownQuantityName )
      call output ( " quantity was not found in the vector", advance='yes' )
    case ( source_not_in_db )
      call output ( " source was not found in the db.", advance='yes' )
    case ( zeroProfilesFound )
      call output ( " command found zero profiles.", advance='yes' )
    case ( zeroGeodSpan )
      call output ( " command found zero geod. ang. span.", advance='yes' )
    case ( vectorWontMatchL2GP )
      call output ( " command found no match of vetor and L2GP.", advance='yes' )
    case ( cantFillFromL2AUX )
      call output ( " command could not be filled from L2AUX.", advance='yes' )
    case ( cantFillFromL1B )
      call output ( " command could not be filled from L1B.", advance='yes' )
    case ( errorReadingL1B )
      call output ( " L1B file could not be read.", advance='yes' )
    case ( vectorWontMatchPDef )
      call output ( " command found new and prev. vectors unmatched.", advance='yes' )
    case ( numInstancesisZero )
      call output ( " command has zero instances.", advance='yes' )
    case ( numSurfsisZero )
      call output ( " command has zero surfaces.", advance='yes' )
    case ( objIsFullRank3 )
      call output ( " command array is full rank 3.", advance='yes' )
    case ( otherErrorInFillVector )
      call output ( " command caused an error in FillVector.", advance='yes' )
    case ( n1_is_zero )
      call output ( " command caused an n1=0 error in squeeze.", advance='yes' )
    case ( n2_is_zero )
      call output ( " command caused an n2=0 error in squeeze.", advance='yes' )
    case ( n3_is_zero )
      call output ( " command caused an n3=0 error in squeeze.", advance='yes' )
    case ( m1_too_small )
      call output ( " command caused a m1 too small error in squeeze.", advance='yes' )
    case ( m2_too_small )
      call output ( " command caused a m2 too small error in squeeze.", advance='yes' )
    case ( not_permutation )
      call output ( " command caused an illegal permutation in squeeze.", advance='yes' )
    case ( allocation_err )
      call output ( " command caused an allocation error in squeeze.", advance='yes' )
    case ( deallocation_err )
      call output ( " command caused an deallocation error in squeeze.", advance='yes' )
    case ( noExplicitValuesGiven )
      call output ( " no explicit values given for explicit fill.", advance='yes' )
    case ( noSourceGridGiven )
      call output ( " no sourceGrid field given for gridded fill.", advance='yes' )
    case ( noSourceL2GPGiven )
      call output ( " no sourceL2GP field given for L2GP fill.", advance='yes' )
    case ( noSourceL2AUXGiven )
      call output ( " no sourceL2AUX field given for L2AUX fill.", advance='yes' )
    case ( noSourceQuantityGiven )
      call output ( " no sourceQuantity field given for vector fill.", advance='yes' )
    case ( invalidExplicitFill )
      call output ( " has inappropriate dimensionality for explicit fill.", advance='yes' )
    case ( needTempREFGPH )
      call output ( " needs temperatureQuantity and refGPHquantity.", advance='yes' )
    case ( needH2O )
      call output ( " needs H2OQuantity.", advance='yes' )
    case ( needGeocAltitude )
      call output ( " needs geocAltitudeQuantity.", advance='yes' )
    case ( badTemperatureQuantity )
      call output ( " temperatureQuantity is not temperature", advance='yes' )
    case ( badREFGPHQuantity )
      call output ( " refGPHQuantity is not refGPH", advance='yes' )
    case ( badGeocAltitudeQuantity )
      call output ( " geocAltitudeQuantity is not geocAltitude", advance='yes' )
    case ( nonConformingHydrostatic )
      call output ( " quantities needed for hydrostatic fill do not conform", advance='yes' )
    case ( badUnitsForMaxIterations )
      call output ( " maxIterations should be dimensionless", advance='yes' )
    case ( noSpecialFill )
      call output ( " invalid special fill", advance='yes' )
    case ( badlosvelfill )
      call output ( " incomplete/incorrect information for los velocity", advance='yes' )
    case default
      call output ( " command caused an unrecognized programming error", advance='yes' )
    end select
    if( present(ExtraMessage)) then
      call output(ExtraMessage)
    endif
  end subroutine ANNOUNCE_ERROR


end module Fill
!=============================================================================

!
! $Log$
! Revision 2.35  2001/04/07 00:12:05  pwagner
! Calls trace_end if needed before every return
!
! Revision 2.34  2001/04/05 23:45:39  pwagner
! Deleted all MLSMessages
!
! Revision 2.33  2001/03/29 19:12:40  livesey
! Added gridded data fill
!
! Revision 2.32  2001/03/15 23:28:23  livesey
! Bug fix
!
! Revision 2.31  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.30  2001/03/15 21:12:11  livesey
! Added special fill for losVel, and dealt with new MLSSignals_m
!
! Revision 2.29  2001/03/15 18:40:38  livesey
! Added some more l1b fill options.
!
! Revision 2.28  2001/03/14 05:33:39  livesey
! Added snoop option
!
! Revision 2.27  2001/03/07 22:42:23  livesey
! Got pressure guesser to work
!
! Revision 2.26  2001/03/06 22:41:07  livesey
! New L2AUX stuff
!
! Revision 2.25  2001/03/06 00:34:46  livesey
! Regular commit.
!
! Revision 2.24  2001/03/05 01:20:14  livesey
! Regular commit, hydrostatic stuff in place.
!
! Revision 2.23  2001/03/03 05:54:29  livesey
! Started hydrostic stuff
!
! Revision 2.22  2001/03/03 00:10:14  livesey
! Removed debuging dump.
!
! Revision 2.21  2001/03/03 00:07:40  livesey
! Added fill from l1b
!
! Revision 2.20  2001/02/27 17:39:03  livesey
! Tidied stuff up a bit.
!
! Revision 2.19  2001/02/27 01:25:15  livesey
! Used new ValidateVectorQuantity routine
!
! Revision 2.18  2001/02/27 00:50:53  livesey
! Made sure verticalCoordinate=l_zeta worked for filling from L2GP
!
! Revision 2.17  2001/02/23 18:16:26  livesey
! Regular commit
!
! Revision 2.16  2001/02/21 01:07:34  livesey
! Got the explicit fill to work.
!
! Revision 2.15  2001/02/08 01:17:41  vsnyder
! Simplify access to abstract syntax tree.
!
! Revision 2.14  2001/01/26 00:11:12  pwagner
! Can fill from prev. defd. vector
!
! Revision 2.13  2001/01/24 23:31:00  pwagner
! Using announce_error, simplified
!
! Revision 2.12  2001/01/10 21:47:45  pwagner
! Chunk bounds determined by geodet.ang.
!
! Revision 2.11  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.10  2001/01/03 17:49:49  pwagner
! Accounts for chunking when filling from old L2GP
!
! Revision 2.9  2000/12/07 00:41:46  pwagner
! added whatquantitynumber
!
! Revision 2.8  2000/12/06 00:01:20  pwagner
! Completed FillOL2AUXData; changed squeeze, nearby
!
! Revision 2.7  2000/12/05 00:40:50  pwagner
! Added FillOL2AUXVector
!
! Revision 2.6  2000/11/30 00:22:52  pwagner
! functions properly moved to read a priori
!
! Revision 2.5  2000/11/16 02:15:25  vsnyder
! Implement timing.
!
! Revision 2.4  2000/11/13 23:02:21  pwagner
! Adapted for rank2 vectorsModule
!
! Revision 2.3  2000/10/06 22:18:47  pwagner
! Fills from old l2gp data
!
! Revision 2.2  2000/09/11 19:52:51  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!
!
