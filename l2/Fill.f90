! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
!=============================================================================

  USE L2GPData
  use GriddedData, only: GriddedData_T
  use INIT_TABLES_MODULE, only: F_SOURCE, S_TIME, S_VECTOR ! Later, s_Fill
                                                           ! will be added
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, NameLen
  use MLSStrings, only: lowercase
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
    & qtyTemplates , L2GPDatabase)

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
    INTEGER :: VectorDBSize
    CHARACTER (LEN=NameLen) :: vectorNameString, templateNameString
    CHARACTER (LEN=NameLen) :: sourceNameString, quantityNameString


    ! These should be moved into open_init some day soon
!    INTEGER :: mlspcf_ol2gp_start=22000
!    INTEGER :: mlspcf_ol2gp_end=22050

!    INTEGER :: OL2FileHandle
    INTEGER :: qtiesStart
    double precision :: T1, T2     ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )

    ! Logical id of file(s) holding old L2GP data
!    OL2FileHandle = mlspcf_ol2gp_start

    ! starting quantities number for *this* vector; what if we have more?
    qtiesStart = 1


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
          ! We can only fill vectors from an old l2gp file for now
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
          CALL FillOL2GPVector(L2GPDatabase(1), qtyTemplates, &
               & vectors(vectorIndex), quantityNameString, qtiesStart)

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
      call output ( t2 - t1, advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Fill

! =====     Private Procedures     =====================================

!=============================== FillVector ==========================
SUBROUTINE FillVector(inPointer, Vector, pointerType, vectorType, NumQtys, &
     & numChans, numSurfs, numInstances, qtiesStart)
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
REAL(r8), POINTER, DIMENSION(:,:,:) ::  inPointer
TYPE(Vector_T), INTENT(OUT) ::          Vector
INTEGER, INTENT(IN) ::                  NumQtys, qtiesStart
INTEGER, INTENT(IN) ::                  numInstances, numSurfs, numChans

! Private
INTEGER ::                              qty, IERR
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
	CALL squeeze(inPointer, Vector%quantities(qty)%values, IERR)
   ENDDO
CASE default
   ! FillVector not yet written to handle these cases
END SELECT

END SUBROUTINE FillVector

!=============================== FillOL2GPVector ==========================
SUBROUTINE FillOL2GPVector(OldL2GPData, qtyTemplates, &
& Output, QuantityName, qtiesStart)
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

! Local variables
!::::::::::::::::::::::::: LOCALS :::::::::::::::::::::
!TYPE(L2GPData_T) ::                               L2GPData
! INTEGER ::                                        Qty
type (QuantityTemplate_T) ::                      OQTemplate	! Output Quantity Template
INTEGER ::                                        i
INTEGER ::                                        ONTimes
INTEGER ::                                        alloc_err
INTEGER ::                                        noL2GPValues=1
LOGICAL ::                                        TheyMatch

OQTemplate = qtyTemplates(Output%TEMPLATE%QUANTITIES(qtiesStart))
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
TheyMatch = OQTemplate%noInstances .EQ. OldL2GPData%nTimes
TheyMatch = TheyMatch .AND. OQTemplate%stacked
IF(TheyMatch) THEN
   DO i = 1, OldL2GPData%nTimes
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
   CALL FillVector(OldL2GPData%l2gpValue, Output, &
        & 'l2gp', 'l2gp', NoL2GPValues, &
        & OldL2GPData%nFreqs, &
        & OldL2GPData%nLevels, &
        & OldL2GPData%nTimes, &
        & qtiesStart)
ELSE
   CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Old and new L2GPData don''t match in times or geolocations')
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

!=============================== nearby ==========================
FUNCTION nearby(x, y)
!=============================== nearby ==========================
! This functions returns TRUE if x and y are "nearby" compared
! with their moduli (proportionally small) and an overall tolerance
REAL(r8), INTENT(IN) ::              x, y

! result
LOGICAL ::                           nearby

! Private
REAL(r8) ::                          small=1.D-32
! REAL(r8) ::                          tolerance=1.D-32

IF(x*y /= 0) THEN
   nearby = ABS(x-y) .LT. small*SQRT(ABS(x))*SQRT(ABS(y))
ELSEIF(x /= 0) THEN
   nearby = ABS(x-y) .LT. small*ABS(x)
ELSE
   nearby=.TRUE.
ENDIF

END FUNCTION nearby

!=============================== squeeze ==========================
SUBROUTINE squeeze(source, sink, IERR)
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
    ! (A future improvement might take as optional arguments
    !  integer arrays source_shape, sink_shape, 
    !  or else shape-params n1, n2, m1)
    !--------Argument--------!
    REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: source
    REAL(r8), DIMENSION(:,:), INTENT(OUT)  :: sink
    INTEGER, INTENT(OUT)                   :: IERR

    !----------Local vars----------!
    ! Error codes
    INTEGER               :: n1_is_zero = 1
    INTEGER               :: n2_is_zero = 2
    INTEGER               :: n3_is_zero = 3
    INTEGER               :: m1_too_small = 4
    INTEGER               :: m2_too_small = 5
    
    INTEGER, DIMENSION(4) :: source_shape
    INTEGER, DIMENSION(4) :: sink_shape
    INTEGER::i,icode,offset
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

    sink = reshape(source, sink_shape(1:2))
    
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
