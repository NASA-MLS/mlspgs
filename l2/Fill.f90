! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
!=============================================================================

  use GriddedData, only: GriddedData_T
  use INIT_TABLES_MODULE, only: F_SOURCE, S_TIME, S_VECTOR ! Later, s_Fill
                                                           ! will be added
  USE L2GPData, only: L2GPData_T
  USE L2AUXData, only: L2AUXData_T, L2AUXDim_None, L2AUXDim_Channel, &
  & L2AUXDim_IntermediateFrequency, L2AUXDim_USBFrequency, L2AUXDim_LSBFrequency, &
  & L2AUXDim_MIF, L2AUXDim_MAF, L2AUXDim_GeodAngle
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, NameLen, LineLen, MLSChunk_T, R8
  use MLSMessageModule, only: MLSMSG_Error, MLSMessage
!  use MLSStrings, only: lowercase
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use string_table, only: get_string
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, &
    & SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_DOT
  use VectorsModule, only: AddVectorToDatabase, CreateVector, Dump, Vector_T, &
    & VectorTemplate_T

  implicit none
  private
  public :: MLSL2Fill
  
  ! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

  ! Error codes for "announce_error"
  integer, parameter :: WRONG_NUMBER = 1     ! of fields of a VECTOR command

  integer, parameter :: s_Fill = 0   ! to be replaced by entry in init_tables_module
  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$id: fill.f90,v 1.1 2000/01/21 21:04:06 livesey Exp $"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module performs the Fill operation in the Level 2 software.  
  ! This takes a vector template, and creates and fills an appropriate vector

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( root, l1bInfo, aprioriData, vectorTemplates, vectors, &
    & qtyTemplates , L2GPDatabase, L2AUXDatabase, chunks, chunkNo )

  ! This is the main routine for the module.  It parses the relevant lines
  ! of the l2cf and works out what to do.

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (L1BInfo_T), intent(in) :: l1bInfo
    type (GriddedData_T), dimension(:), pointer :: aprioriData
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (Vector_T), dimension(:), pointer :: vectors
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: L2GPDatabase
    TYPE (L2AUXData_T), DIMENSION(:), POINTER :: L2AUXDatabase
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    integer, intent(in) :: chunkNo

    ! Local variables
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Definitely n_named
    type (Vector_T) :: newVector
    integer :: SON                 ! Of root, an n_spec_args or a n_named
    integer :: templateIndex       ! In the template database
    integer :: vectorIndex         ! In the vector database
    integer :: vectorName          ! Sub-rosa index
    integer :: quantityName        ! Sub-rosa index
    integer :: sourceName          ! Sub-rosa index
    integer :: l2Index             ! Where source is among l2gp or l2aux database
    INTEGER :: VectorDBSize
    CHARACTER (LEN=NameLen) :: vectorNameString, templateNameString
    CHARACTER (LEN=NameLen) :: sourceNameString, quantityNameString
    LOGICAL :: is_l2gp, is_l2aux

    ! These should be moved into open_init some day soon
!    INTEGER :: mlspcf_ol2gp_start=22000
!    INTEGER :: mlspcf_ol2gp_end=22050

!    INTEGER :: OL2FileHandle
    INTEGER :: qtiesStart
    CHARACTER (LEN=LineLen) ::              msr
    REAL :: T1, T2     ! for timing
    logical :: TIMING

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

      ! Node_id(key) is now n_spec_args.

      select case( decoration(subtree(1,decoration(subtree(1,key)))) )
      case ( s_vector )
        if ( nsons(key) /= 2 ) call announce_error ( son, wrong_number )
        templateIndex = decoration(decoration(subtree(2,subtree(2,key))))

        ! Create the vector, and add it to the database.

        call decorate ( key, AddVectorToDatabase ( vectors, &
          & CreateVector ( vectorName, vectorTemplates(templateIndex), &
            & qtyTemplates ) ) )

        ! That's the end of the create operation

     case ( s_fill )
          ! We can for now only fill vectors from either:
        ! old l2gp file
        ! old l2aux file
        !          vectorName=""
        !          sourceName=""
        !          quantityName=""

          ! A Fill instruction, this contains just a vector name and quantity 
          ! and a source name assuming the same quantity
          ! e.g., Fill, state.quantity, source=oldState.quantity

        if ( nsons(key) < 3 ) call announce_error ( son, wrong_number )
        IF(node_id(subtree(2, key)) /= n_dot) then
           ! In case dotless
           vectorIndex = decoration(decoration(subtree(2,subtree(2,key))))
           quantityName=0
        ELSE
           ! Assume form : x.y
           J = subtree(2, subtree(2, key))
           vectorName=sub_rosa(subtree(2, J))     ! x
           quantityName=sub_rosa(subtree(3, J))   ! y
           vectorIndex = decoration(subtree(2, J))
        ENDIF

        DO J=3, nsons(key)
           SELECT CASE( decoration(subtree(1,decoration(subtree(J,key)))) )
              CASE(f_source)
                 ! source=vector.quantity
                 sourceName=sub_rosa(subtree(J, subtree(2, key)))
                 IF(quantityName==0) THEN
                    quantityName=sub_rosa(subtree(3, subtree(J, key)))   ! y
                 ENDIF
!              CASE(f_method)
!              CASE(f_other)
              CASE DEFAULT ! Can't get here if tree_checker worked correctly
           END SELECT
        ENDDO

        CALL get_string(vectorName, vectorNameString)
        CALL get_string(quantityName, quantityNameString)
        CALL get_string(sourceName, sourceNameString)

!          vectorIndex=LinearSearchStringArray(vectors%name,&
!               & vectorName,caseInsensitive=.FALSE.)

         ! compute what quantity number the vector quantity corresponds to
         qtiesStart = whatQuantityNumber(vectorName, quantityName, &
         & vectors(vectorIndex)%Template)

         IF(qtiesStart < 0) THEN
              msr = 'Quantity Name ' // quantityNameString // &
              & ' not found among quantities assoc. with ' // vectorNameString
              CALL MLSMessage(MLSMSG_Error, ModuleName, &
                     & msr)
         ENDIF
! Is our source of type l2gp or l2aux?
         is_l2gp = .FALSE.
         is_l2aux = .FALSE.
 
         l2Index = 1
         DO WHILE (.NOT. is_l2gp .AND. l2Index.LE.SIZE(L2GPDatabase))
             IF(sourceName.EQ.L2GPDatabase(l2Index)%nameIndex) THEN
                 is_l2gp = .TRUE.
                 exit
             ENDIF
             l2Index = l2Index + 1
         ENDDO
 
         IF(.NOT. is_l2gp) THEN
                l2Index = 1
                DO WHILE (.NOT. is_l2aux .AND. l2Index.LE.SIZE(L2AUXDatabase))
                    IF(sourceName.EQ.L2AUXDatabase(l2Index)%Name) THEN
                        is_l2aux = .TRUE.
                        exit
                    ENDIF
                    l2Index = l2Index + 1
                ENDDO
         ENDIF

! Fill
         IF(is_l2gp) THEN
                CALL FillOL2GPVector(L2GPDatabase(l2Index), qtyTemplates, &
                      & vectors(vectorIndex), quantityNameString, qtiesStart, &
                      & chunkNo)

         ELSEIF(is_l2aux) THEN
                CALL FillOL2AUXVector(L2AUXDatabase(l2Index), qtyTemplates, &
                      & vectors(vectorIndex), quantityNameString, qtiesStart, &
                      & chunkNo)

         ELSE
                CALL MLSMessage(MLSMSG_Error, ModuleName, &
                     & 'source file in neither l2gp nor l2aux databases')
         ENDIF
          ! This is *not* going to work if you have more than one vector

          qtiesStart = qtiesStart+1

      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

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
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Fill

! =====     Private Procedures     =====================================

!=============================== FillVector ==========================
SUBROUTINE FillVector(inPointer, Vector, pointerType, vectorType, NumQtys, &
     & numChans, numSurfs, numInstances, qtiesStart, dim_order)
!=============================== FillVector ==========================

! Fill the vector Vector with values taken from the array pointed to by inPointer
! in a manner that depends on their respective types:
!
!        pointerType        vectorType            operation
!          l2gp              l2gp       vector(:,:,:) = pointer(:,:,:)

! It is assumed that the rank3 Pointer is filled from
! Pointer(1, 1, 1) to Pointer(numChans, numSurfs, numInstances)
!
! With the redefinition of Vector%quantities to be a rank2 object
! we are faced with a problem: how to fill a rank 2 object from a rank 3 one?
! If the rank 3 object is full, then it use the following trick:
! Vector(1:numChans*numSurfs, 1:numInstances) = 
!                   Pointer(1:numChans, 1:numSurfs, 1:numInstances)

CHARACTER (Len=*), INTENT(IN) ::        pointerType, vectorType
! REAL(r8), POINTER, DIMENSION(:,:,:) ::  inPointer
REAL(r8), DIMENSION(:,:,:) ::           inPointer
TYPE(Vector_T), INTENT(OUT) ::          Vector
INTEGER, INTENT(IN) ::                  NumQtys, qtiesStart
INTEGER, INTENT(IN) ::                  numInstances, numSurfs, numChans
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: dim_order

! Private
CHARACTER (LEN=LineLen) ::              msr
INTEGER ::                              qty, IERR
CHARACTER (LEN=4) ::                    qtyChar, IERRChar
! Error codes
INTEGER ::                              numInstancesisZero=1
INTEGER ::                              numSurfsisZero=2
INTEGER ::                              numChansisZero=3
INTEGER ::                              objIsFullRank3=4
INTEGER ::                              ErrorInFillVector=-999

! Sanity checks:
IF(numChans.EQ.0) THEN
	call announce_error ( ErrorInFillVector, numChansisZero )
        RETURN
ELSEIF(numSurfs.EQ.0) THEN
	call announce_error ( ErrorInFillVector, numSurfsisZero )
        RETURN
ELSEIF(numInstances.EQ.0) THEN
	call announce_error ( ErrorInFillVector, numInstancesisZero )
        RETURN
ENDIF

Vector%template%NoQuantities = NumQtys
SELECT CASE (PointerType(:4) // VectorType(:4))
CASE ('l2gpl2gp')
   DO qty = qtiesStart, qtiesStart - 1 + NumQtys
     WRITE(qtyChar, '(I4)') qty
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
	CALL squeeze(IERR, inPointer, Vector%quantities(qty)%values)
     WRITE(IERRChar, '(I4)') IERR
     msr = 'Error #' // IERRChar // ' in squeezing l2gp into quantity ' // qtyChar
     IF(IERR /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & msr)
   ENDDO
CASE ('l2aul2au')
   DO qty = qtiesStart, qtiesStart - 1 + NumQtys
     WRITE(qtyChar, '(I4)') qty
     Vector%quantities(qty)%template%noChans = numChans
     Vector%quantities(qty)%template%noSurfs = numSurfs
     Vector%quantities(qty)%template%noInstances = numInstances
     IF(PRESENT(dim_order)) THEN
        CALL squeeze(IERR, inPointer, Vector%quantities(qty)%values, dim_order)
     ELSE    
        CALL squeeze(IERR, inPointer, Vector%quantities(qty)%values)
     ENDIF
     WRITE(IERRChar, '(I4)') IERR
     msr = 'Error #' // IERRChar // ' in squeezing l2aux into quantity ' // qtyChar
     IF(IERR /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & msr)
   ENDDO
CASE default
   ! FillVector not yet written to handle these cases
END SELECT

END SUBROUTINE FillVector

!=============================== FillOL2GPVector ==========================
SUBROUTINE FillOL2GPVector(OldL2GPData, qtyTemplates, &
& Output, QuantityName, qtiesStart, chunkNo)
!=============================== FillOL2GPVector ==========================

! If the times, pressures, and geolocations match,
! fill the vector Output with values taken from the appropriate quantity in
! Old L2GP vector OldL2GPData
! 

INTEGER, INTENT(IN) ::                            qtiesStart
TYPE(L2GPData_T) ::                               OldL2GPData
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
!TYPE(L2GPData_T), DIMENSION(:), POINTER ::        L2GPDatabase
TYPE(Vector_T), INTENT(INOUT) ::                  Output
CHARACTER*(*), INTENT(IN) ::                      QuantityName
! If done chunk-by-chunk, the following is the chunk number
! Otherwise, chunkNo should = -1
    integer, intent(in) :: chunkNo

! Local variables
!::::::::::::::::::::::::: LOCALS :::::::::::::::::::::
!TYPE(L2GPData_T) ::                               L2GPData
! INTEGER ::                                        Qty
type (QuantityTemplate_T) ::                      OQTemplate	! Output Quantity Template
INTEGER ::                                        i
INTEGER ::                                        ONTimes
INTEGER ::                                        alloc_err
INTEGER ::                                        firstProfile, lastProfile
INTEGER ::                                        noL2GPValues=1
LOGICAL ::                                        TheyMatch

OQTemplate = qtyTemplates(Output%TEMPLATE%QUANTITIES(qtiesStart))
! Chunk-by-chunk, or all chunks at once?
IF(chunkNo == -1) THEN
	firstProfile = 1
   lastProfile = OldL2GPData%nTimes
ELSE
	firstProfile = 1
   lastProfile = OldL2GPData%nTimes
	DO i=1, OldL2GPData%nTimes
        	IF(chunkNo == OldL2GPData%chunkNumber(i)) THEN
				firstProfile = MIN(firstProfile, i)
				lastProfile = MAX(lastProfile, i)
         ENDIF
   ENDDO
ENDIF
!
! Read the old L2GP file for QuantityName
!CALL ReadL2GPData(OL2FileHandle, TRIM(LowerCase(QuantityName)), &
!     & OldL2GPData, ONTimes)
! Allocate space for current l2gpdata
!ALLOCATE(l2gpData%pressures(OldL2GPData%nLevels),&
!     & l2gpData%latitude(ONTimes), &
!     & l2gpData%longitude(ONTimes), & 
!     & l2gpData%time(ONTimes), &
!     & l2gpData%solarTime(ONTimes), &
!     & l2gpData%solarZenith(ONTimes), &
!     & l2gpData%losAngle(ONTimes), &
!     & l2gpData%geodAngle(ONTimes), &
!     & l2gpData%chunkNumber(ONTimes), &
!     & l2gpData%frequency(OldL2GPData%nFreqs), &
!     & l2gpData%l2gpValue(OldL2GPData%nFreqs, OldL2GPData%nLevels, ONTimes), &
!     & l2gpData%l2gpPrecision(OldL2GPData%nFreqs, OldL2GPData%nLevels, ONTimes), &
!     & l2gpData%status(ONTimes), l2gpData%quality(ONTimes), &
!     & STAT=alloc_err)
!IF(alloc_err /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
!     & 'Failed to allocate temp l2gpData')
! Check that times, etc. match
!L2GPData    = L2GPDatabase(1)
TheyMatch = OQTemplate%noInstances .EQ. (lastProfile-firstProfile+1)
TheyMatch = TheyMatch .AND. OQTemplate%stacked
IF(TheyMatch) THEN
   DO i = firstProfile, lastProfile
      IF( &
        &   nearBy(OldL2GPData%latitude(i), OQTemplate%geodLat(1,i)) &
        & .AND. &
        &   nearBy(OldL2GPData%longitude(i), OQTemplate%lon(1,i)) &
        & .AND. &
        &   nearBy(OldL2GPData%solarTime(i), OQTemplate%solarTime(1,i)) &
        &) THEN
       ELSE
          TheyMatch = .FALSE.
          EXIT
      ENDIF
   END DO
ENDIF
IF(TheyMatch .AND. OldL2GPData%NLevels.EQ.OQTemplate%noSurfs) THEN
   DO i = 1, OldL2GPData%NLevels
      IF( &
        &   nearBy(OldL2GPData%pressures(i), OQTemplate%surfs(i,1)) &
        & ) THEN
     ELSE
       TheyMatch = .FALSE.
       EXIT
    ENDIF
   END DO
ENDIF
IF(TheyMatch) THEN
!   DO Qty=1, 
!   CALL put_3d_view(Output, &
!        & OldL2GPData%quantities(qty)%template%noFreqs, &
!        & OldL2GPData%quantities(qty)%template%noSurfs, &
!        & OldL2GPData%quantities(qty)%template%noInstances, &
!        & OldL2GPData%quantities(qty)%values)
   CALL FillVector(OldL2GPData%l2gpValue(:, :, firstProfile:lastProfile), &
        & Output, &
        & 'l2gp', 'l2gp', NoL2GPValues, &
        & OldL2GPData%nFreqs, &
        & OldL2GPData%nLevels, &
        & lastProfile-firstProfile+1, &
        & qtiesStart)
ELSE
   CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Vector and old L2GP do not match in times or geolocations')
ENDIF
!DEALLOCATE(l2gpData%pressures, &
!     & l2gpData%latitude, l2gpData%longitude, l2gpData%time, &
!     & l2gpData%solarTime, l2gpData%solarZenith, l2gpData%losAngle, &
!     & l2gpData%geodAngle, l2gpData%chunkNumber, &
!     & l2gpData%l2gpValue, l2gpData%frequency, &
!     & l2gpData%l2gpPrecision, l2gpData%status, l2gpData%quality, &
!     & STAT=alloc_err)
!IF(alloc_err /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
!     & 'Failed to deallocate temp l2gpData')
!
!DEALLOCATE(OldL2GPData%pressures, &
!     & OldL2GPData%latitude, l2gpData%longitude, OldL2GPData%time, &
!     & OldL2GPData%solarTime, OldL2GPData%solarZenith, OldL2GPData%losAngle, &
!     & OldL2GPData%geodAngle, OldL2GPData%chunkNumber, &
!     & OldL2GPData%l2gpValue, OldL2GPData%frequency, &
!     & OldL2GPData%l2gpPrecision, OldL2GPData%status, OldL2GPData%quality,&
!     & STAT=alloc_err)
!IF(alloc_err /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, &
!     & 'Failed to deallocate temp OldL2GPData')
END SUBROUTINE FillOL2GPVector

!=============================== FillOL2AUXVector ==========================
SUBROUTINE FillOL2AUXVector(OldL2AUXData, qtyTemplates, &
& Output, QuantityName, qtiesStart, chunkNo)
!=============================== FillOL2AUXVector ==========================

! If the times, pressures, and geolocations match,
! fill the vector Output with values taken from the
! Old L2AUX vector OldL2AUXData
! 

INTEGER, INTENT(IN) ::                            qtiesStart
TYPE(L2AUXData_T) ::                               OldL2AUXData
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates
!TYPE(L2GPData_T), DIMENSION(:), POINTER ::        L2GPDatabase
TYPE(Vector_T), INTENT(INOUT) ::                  Output
CHARACTER*(*), INTENT(IN) ::                      QuantityName
    integer, intent(in) :: chunkNo

! Local variables
!::::::::::::::::::::::::: LOCALS :::::::::::::::::::::
type (QuantityTemplate_T) ::                      OQTemplate	! Output Quantity Template
INTEGER ::                                        OQType	! Ouptut quantity type
INTEGER ::                                        OQNVals	! No. Ouptut quantity vals
INTEGER ::                                        ChanDim	! which dim no. is channel
INTEGER ::                                        IntFreqDim	! which dim no. is int. frq.
INTEGER ::                                        USBDim	! which dim no. is USB
INTEGER ::                                        LSBDim	! which dim no. is LSB
INTEGER ::                                        MAFDim	! which dim no. is MAF
INTEGER ::                                        MIFDim	! which dim no. is MIF
INTEGER ::                                        L2NVals	! No. l2 aux vals
INTEGER ::                                        i
LOGICAL ::                                        TheyMatch
! This will let us re-order the vector dimensions differently from the l2aux
INTEGER, DIMENSION(3) ::                          dim_order

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
DO i=1, OldL2AUXData%noDimensionsUsed
    L2NVals = L2NVals*OldL2AUXData%dimensions(i)%noValues
    SELECT CASE(OldL2AUXData%dimensions(i)%dimensionFamily)
    CASE(L2AUXDim_Channel)
        ChanDim = i
        dim_order(1) = i
    CASE(L2AUXDim_IntermediateFrequency)
        IntFreqDim = i
    CASE(L2AUXDim_USBFrequency)
        USBDim = i
    CASE(L2AUXDim_LSBFrequency)
        LSBDim = i
    CASE(L2AUXDim_MIF)
        MIFDim = i
    CASE(L2AUXDim_MAF)
        MAFDim = i
    CASE DEFAULT	! We are not yet interested in these dimensions
    END SELECT
    
ENDDO

! Check that the dimensions match up and are of the same family
TheyMatch = OQNVals == L2NVals
IF(OQTemplate%minorFrame) THEN
	TheyMatch = TheyMatch .AND. MIFDim > 0
ENDIF
IF(OQTemplate%noChans > 0) THEN
	TheyMatch = TheyMatch .AND. &
        & (  ChanDim > 0   )
ENDIF

IF(TheyMatch) THEN
   CALL FillVector(OldL2AUXData%values, Output, &
        & 'l2aux', 'l2aux', 1, &
        & MAX(OldL2AUXData%dimensions(1)%noValues, 1), &
        & MAX(OldL2AUXData%dimensions(2)%noValues, 1), &
        & MAX(OldL2AUXData%dimensions(3)%noValues, 1), &
        & qtiesStart)
ELSE
   CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Vector and old L2AUX do not match')
ENDIF

END SUBROUTINE FillOL2AUXVector

!=============================== nearby ==========================
FUNCTION nearby(x, y, inRelSmall, inOverallTol)
!=============================== nearby ==========================
! This functions returns TRUE if x and y are "nearby" comparing either
! their relative difference (proportionally small) 
! or their differences with an overall tolerance
REAL(r8), INTENT(IN) ::              x, y
REAL(r8), INTENT(IN), OPTIONAL ::    inRelSmall
REAL(r8), INTENT(IN), OPTIONAL ::    inOverallTol

! result
LOGICAL ::                           nearby

! Private
REAL(r8) ::                          small=1.D-32
REAL(r8) ::                          tolerance=0.D0

IF(PRESENT(inRelSmall)) THEN
    small = inRelSmall
ENDIF

IF(PRESENT(inOverallTol)) THEN
    tolerance = inOverallTol
ENDIF

IF(ABS(x-y) <= tolerance) THEN
   nearby = .TRUE.
ELSEIF(x*y /= 0) THEN
   nearby = ABS(x-y) .LT. small*SQRT(ABS(x))*SQRT(ABS(y))
ELSEIF(x /= 0) THEN
   nearby = ABS(x-y) .LT. small*ABS(x)
ELSE
   nearby=.TRUE.
ENDIF

END FUNCTION nearby

!=============================== whatQuantityNumber ==========================
FUNCTION whatQuantityNumber(x, y, xVectorTemplate)
!=============================== whatQuantityNumber ==========================
! This functions returns the quantity number of quantity y in vector x
! where x and y are sub-rosa indexes of their actual names;
! i.e., they are specified in the cf as name[x].name[y]
! where name[x] = get_char(x), name[y] = get_char[y]
!
! If the quantity is not found, it returns FAILED (-1)
INTEGER, INTENT(IN) ::               x, y
type (VectorTemplate_T)  ::          xVectorTemplate

! result
INTEGER ::                           whatQuantityNumber

! Private
INTEGER, PARAMETER ::                FAILED=-1
INTEGER              ::              qty

whatQuantityNumber = FAILED
DO qty = 1, xVectorTemplate%NoQuantities
    IF(xVectorTemplate%Name == y) THEN
        whatQuantityNumber = qty
        EXIT
    ENDIF
ENDDO

END FUNCTION whatQuantityNumber

!=============================== squeeze ==========================
SUBROUTINE squeeze(IERR, source, sink, source_order)
!=============================== squeeze ==========================
    ! takes a rank 3 object source and returns a rank2 object sink
    ! source(1..n1, 1..n2, 1..n3) -> sink(1..n1*n2, 1..n3)
    ! unless it can't--then it returns iERR /= 0
    ! One reason it may fail: shape of sink too small
    !
    ! Assuming that shape(source) = {n1, n2, n3}
    !     =>          shape(sink) = {m1, m2}
    ! then we must further assume (else set IERR)
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
    REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: source
    REAL(r8), DIMENSION(:,:), INTENT(OUT)  :: sink
    INTEGER, INTENT(OUT)                   :: IERR
    INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: source_order

    !----------Local vars----------!
    ! Error codes
    INTEGER               :: n1_is_zero = 1
    INTEGER               :: n2_is_zero = 2
    INTEGER               :: n3_is_zero = 3
    INTEGER               :: m1_too_small = 4
    INTEGER               :: m2_too_small = 5
    INTEGER               :: not_permutation = 6
    INTEGER               :: allocation_err = 7
    INTEGER               :: deallocation_err = 8
    
    INTEGER, DIMENSION(4) :: source_shape
    INTEGER, DIMENSION(4) :: sink_shape
    INTEGER, DIMENSION(4) :: temp_shape
    INTEGER::i,icode,offset
    REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: temp
    !----------Executable part----------!
    source_shape(1:3) = shape(source)
    sink_shape(1:2) = shape(sink)

    IF (source_shape(1) == 0) THEN
    	IERR = n1_is_zero
        RETURN
    ELSEIF (source_shape(2) == 0) THEN
    	IERR = n2_is_zero
        RETURN
    ELSEIF (source_shape(3) == 0) THEN
    	IERR = n3_is_zero
        RETURN
    ELSEIF (sink_shape(1) < source_shape(1)*source_shape(2)) THEN
    	IERR = m1_too_small
        RETURN
    ELSEIF (sink_shape(2) < source_shape(3)) THEN
    	IERR = m2_too_small
        RETURN
    ELSE
    	IERR = 0
    ENDIF

    IF(PRESENT(source_order)) THEN
!        Check that source_order is a legal permutation of (1 2 3)
!        using trick: sum and product must each equal 6
        IF(source_order(1)+source_order(2)+source_order(3) /= 6 &
        & .OR. &
        & source_order(1)*source_order(2)*source_order(3) /= 6 ) THEN
             IERR=6
             RETURN
        ENDIF
        DO i=1, 3
           temp_shape(i) = source_shape(source_order(i))
        ENDDO
        ALLOCATE(temp(temp_shape(1), temp_shape(2), temp_shape(3)), &
        & STAT=IERR)
        IF(IERR /= 0) THEN
        	IERR=allocation_err
                RETURN
        ENDIF
        temp = reshape(source, temp_shape(1:3), order=source_order(1:3))
        sink = reshape(temp, sink_shape(1:2))
        DEALLOCATE(temp, &
        & STAT=IERR)
        IF(IERR /= 0) THEN
        	IERR=deallocation_err
                RETURN
        ENDIF
    ELSE
        sink = reshape(source, sink_shape(1:2))
    ENDIF
    
END SUBROUTINE squeeze

!
!
  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( wrong_number )
      call output ( "The " );
      call dump_tree_node ( where, 0 )
      call output ( " command does not have exactly one field.", advance='yes' )
    end select
  end subroutine ANNOUNCE_ERROR

!=============================================================================
end module Fill
!=============================================================================

!
! $Log$
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
