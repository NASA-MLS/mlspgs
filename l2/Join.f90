! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Join                     ! Join together chunk based data.
!=============================================================================

  ! This module performs the 'join' task in the MLS level 2 software.

  use INIT_TABLES_MODULE, only: F_COMPAREOVERLAPS, F_FILE, F_OUTPUTOVERLAPS, F_SOURCE, &
    & F_UNPACKOUTPUT, FIELD_FIRST, FIELD_INDICES, FIELD_LAST, L_PRESSURE, L_NONE, &
    & L_TRUE, L_ZETA, S_L2AUX, S_L2GP, S_TIME
  use L2AUXData, only: AddL2AUXToDatabase, ExpandL2AUXDataInPlace, &
    & L2AUXData_T, L2AUXDim_Channel, L2AUXDim_geodAngle, &
    & L2AUXDim_IntermediateFrequency, L2AUXDim_LSBFrequency, L2AUXDim_MAF, &
    & L2AUXDim_MIF, L2AUXDim_None, L2AUXDim_USBFrequency, SetupNewL2AUXRecord
  use L2GPData, only: AddL2GPToDatabase, ExpandL2GPDataInPlace, &
    & L2GPData_T, SetupNewL2GPRecord
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: MLSChunk_T, R8
! use MLSL2Common
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use SDPToolkit, only: PGS_S_SUCCESS, PGS_TD_TAItoUTC, PGSTD_E_NO_LEAP_SECS
  use String_Table, only: DISPLAY_STRING
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
    & SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_SET_ONE
  use VectorsModule, only: GetVectorQuantity, GetVectorQtyByTemplateIndex, &
    & ValidateVectorQuantity, Vector_T, VectorValue_T, DUMP
  use Intrinsic, ONLY: L_NONE, L_INSTRUMENTCHANNEL, L_USBFREQUENCY, L_LSBFREQUENCY,&
       L_INTERMEDIATEFREQUENCY

  implicit none
  private
  public :: MLSL2Join

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Parameters for Announce_Error

  integer :: ERROR
  integer, parameter :: NotAllowed=1

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  MLSL2Join  -----

  ! This is the main routine for join.  Most of the time it is fairly simple.
  ! However, for the first time round a little more has to be done as the
  ! routine has to create the l2gp and l2aux structures with the correct size
  ! in order to be able to store all the chunks.

  subroutine MLSL2Join ( root, vectors, l2gpDatabase, l2auxDatabase, chunkNo )

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the JOIN section in the AST
    type (Vector_T), dimension(:), intent(in) :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo

    ! Local variables
    integer :: FIELD               ! Subtree index of "field" node
    integer :: FIELD_INDEX         ! F_..., see Init_Tables_Module
    logical :: GOT_FIELD(field_first:field_last)
    integer :: GSON                ! Son of Key
    integer :: KEY                 ! Index of an L2GP or L2AUX tree
    integer :: KEYNO               ! Index of subtree of KEY
    integer :: mlscfLine
    integer :: NAME                ! Sub-rosa index of name of L2GP or L2AUX
    integer :: SON                 ! A son of ROOT
    integer :: SOURCE              ! Index in AST
    integer :: VALUE               ! Value of a field
    integer :: vectorIndex, quantityIndex
    type (VectorValue_T), pointer :: quantity
    logical :: compareOverlaps, outputOverlaps, unpackOutput
    REAL :: T1, T2     ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Join", root )

    error = 0

    ! We simply loop over the lines in the mlscf

    do mlscfLine = 2, nsons(root)-1     ! Skip name at begin and end of section
      ! Each line represents a different join operation; clear various
      ! flags etc.
      son = subtree(mlscfLine,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        name = sub_rosa(subtree(1,son))
      else
        key = son
        name = 0
      end if

      ! Node_id(key) is now n_spec_args.

      ! ??? Does this need to do anything somewhere ???
      select case( decoration(subtree(1,decoration(subtree(1,key)))) )
      case ( s_l2aux )
      case ( s_l2gp )
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select

      got_field = .false.
      source = null_tree
!     name=mlscfSection%entries(mlscfLine)%mlscfEntryName
      compareOverlaps = .FALSE.
      outputOverlaps = .FALSE.
      unpackOutput = .FALSE.

      ! Loop over the fields of the mlscf line

      do keyNo = 2, nsons(key) ! Skip spec name
        gson = subtree(keyNo,key)
        field = subtree(1,gson)
        if ( node_id(gson) == n_set_one ) then
          value = l_true
        else
          value = decoration(subtree(2,gson))
        end if
        field_index = decoration(field)
        got_field(field_index) = .true.
        select case ( field_index )
        case ( f_source )
          source = subtree(2,gson) ! required to be an n_dot vertex
          vectorIndex = decoration(decoration(subtree(1,source)))
          quantityIndex = decoration(decoration(decoration(subtree(2,source))))
        case ( f_compareoverlaps )
          compareOverlaps = value == l_true
        case ( f_outputoverlaps )
          outputOverlaps = value == l_true
        case ( f_file)
          call announce_error(key,NotAllowed,field_index)
        case ( f_unpackoutput )
          unpackOutput = value == l_true
        case default ! Can't get here if tree_checker worked properly
        end select
      end do

      if ( error > 0 ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Errors in configuration prevent proceeding" )

      quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex),quantityIndex)

      ! Now, depending on the properties of the source we deal with the
      ! vector quantity appropriately.

      if (ValidateVectorQuantity(quantity,coherent=.TRUE.,stacked=.TRUE.,regular=.TRUE.,&
        & verticalCoordinate=(/L_Pressure,L_Zeta,L_None/))) then
        ! Coherent, stacked, regular quantities on pressure surfaces, or
        ! with no vertical coordinate system go in l2gp files.
        call display_string(quantity%template%name)
        call JoinL2GPQuantities ( key, name, quantity, l2gpDatabase, chunkNo )
      else
        ! All others go in l2aux files.
        call JoinL2AUXQuantities ( key, name, quantity, l2auxDatabase, chunkNo, &
          & unpackOutput=unpackOutput )
      endif
    end do

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 ) then
!       call dump ( ??? )
      end if
      call trace_end ( "MLSL2Join" )
    end if
    if ( timing ) call sayTime

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for MLSL2Join =" )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine MLSL2Join

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  Announce_Error  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, FIELDINDEX )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), OPTIONAL :: FIELDINDEX ! Extra information for msg

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
      case ( NotAllowed )
        call output('Field ')
        call display_string(field_indices(fieldIndex))
        call output(' is not allowed in this context',advance='yes')
    end select
  end subroutine ANNOUNCE_ERROR
  ! -----------------------------------------  JoinL2GPQuantities  -----

  ! This routine joins an l2gp line quantity to a database of such quantities.
  ! If this is the first time through, the database is created.

  ! The firstInstance and lastInstance arguments give an optional range of
  ! the instances that we wish to store in the l2gp quantity.  Otherwise, it
  ! defaults to the non overlapped region.

  subroutine JoinL2GPQuantities ( key, name, quantity, l2gpDatabase, &
    & chunkNo, firstInstance, lastInstance )

    ! Dummy arguments
    integer, intent(in) :: KEY     ! spec_args to Decorate with the L2GP index
    integer, intent(in) :: NAME    ! Of the l2gp command
    type (VectorValue_T), intent(in) :: QUANTITY
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    integer, intent(in) :: CHUNKNO
    integer, intent(in), optional :: FIRSTINSTANCE, LASTINSTANCE
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2gp data.

    ! Local variables

    integer :: status,profNo,freqNo
    type (L2GPData_T) :: newL2GP
    type (L2GPData_T), pointer :: thisL2GP
    integer :: index
    integer :: firstProfile,lastProfile ! Profile range in the l2gp to output to
    integer :: noSurfsInL2GP,noFreqsInL2GP
    integer :: useFirstInstance,useLastInstance,noOutputInstances
    logical :: l2gpDataIsNew
    real(r8), dimension(:,:), pointer :: values
    
    if ( toggle(gen) ) call trace_begin ( "JoinL2GPQuantities", key )

    ! If this is the first chunk, we have to setup the l2gp quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.
    l2gpDataIsNew = (.NOT. associated(l2gpDatabase))
    if ( .not. l2gpDataIsNew ) then
      index = decoration(key)
      l2gpDataIsNew = (index>=0)
    end if

    ! Work out what to do with the first and last instance information
    
    if ( PRESENT(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( PRESENT(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    noOutputInstances = useLastInstance-useFirstInstance+1
    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    if ( l2gpDataIsNew ) then
      ! Now create an empty L2GP record with this dimension

      if (any(quantity%template%verticalCoordinate == (/l_Pressure, l_Zeta /) )) then
        noSurfsInL2GP = quantity%template%noSurfs
      else
        noSurfsInL2GP = 0
      end if

      if ( quantity%template%frequencyCoordinate == l_None) then
         noFreqsInL2GP=0
      else
         noFreqsInL2GP=quantity%template%noChans
      endif

      call SetupNewL2GPRecord ( newL2GP, noFreqsInL2GP, noSurfsInL2GP )
      ! Setup the standard stuff, only pressure as it turns out.
      if ( quantity%template%verticalCoordinate == l_Pressure ) &
        & newL2GP%pressures = quantity%template%surfs(:,1)
      if ( quantity%template%verticalCoordinate == l_Zeta ) &
        & newL2GP%pressures = 10.0**(-quantity%template%surfs(:,1))

      ! ??? In later versions we'll need to think about frequency stuff (NJL)

      ! Add it to the database of l2gp quantities
      index = AddL2GPToDatabase ( l2gpDatabase, newL2GP )
      ! Setup the pointer and index to be used later
      call decorate ( key, -index ) ! Remember where it is
      thisL2GP => l2gpDatabase(index)

    else
      ! Setup the index and pointer
      thisL2GP => l2gpDatabase(-index)
    end if

    ! Expand l2gp (initially all zero-size arrays) to take the new information
    call ExpandL2GPDataInPlace ( thisL2GP, &
      & thisL2GP%nTimes+noOutputInstances )

    ! Now copy the information from the quantity to the l2gpData

    ! name is an integer, but L2GP%name is Character data
    thisL2GP%nameIndex = name
    lastProfile=thisL2GP%nTimes
    firstProfile=lastProfile-noOutputInstances+1

    ! Now fill the data, first the geolocation
    thisL2GP%latitude(firstProfile:lastProfile) = &
      & quantity%template%geodLat(1,useFirstInstance:useLastInstance)
    thisL2GP%longitude(firstProfile:lastProfile) = &
      & quantity%template%lon(1,useFirstInstance:useLastInstance)
    thisL2GP%solarTime(firstProfile:lastProfile) = &
      & quantity%template%solarTime(1,useFirstInstance:useLastInstance)
    thisL2GP%solarZenith(firstProfile:lastProfile) = &
      & quantity%template%solarZenith(1,useFirstInstance:useLastInstance)
    thisL2GP%losAngle(firstProfile:lastProfile) = &
      & quantity%template%losAngle(1,useFirstInstance:useLastInstance)
    thisL2GP%geodAngle(firstProfile:lastProfile) = &
      & quantity%template%phi(1,useFirstInstance:useLastInstance)
    thisL2GP%time(firstProfile:lastProfile) = &
      & quantity%template%time(1,useFirstInstance:useLastInstance)
    thisL2GP%chunkNumber(firstProfile:lastProfile)=chunkNo

    ! Now the various data quantities.

    ! For v0.1 we're only thinking about value.  The precision will
    ! come from matrices later in 0.5, and the diagnostics such as status
    ! and quality will come later too (probably 0.5, but maybe 1.0)

    thisL2GP%l2gpValue(:,:,firstProfile:lastProfile)=&
         RESHAPE(quantity%values(:,useFirstInstance:useLastInstance),&
         (/MAX(thisL2GP%nFreqs,1),MAX(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0 ! Later put something here
    thisL2GP%status(firstProfile:lastProfile)='G'
    thisL2GP%quality(firstProfile:lastProfile)=0.0

    if ( toggle(gen) ) call trace_end ( "JoinL2GPQuantities" )
  end subroutine JoinL2GPQuantities

  ! ----------------------------------------  JoinL2AUXQuantities  -----

  ! This subroutine is like the one above, except that the quantities it joins
  ! are destined to go in L2AUX quantities.

  ! The unpackOutput argument warrants some explanation.  If it is set false,
  ! then the data is output with no gaps in the arrays.  If it is true, then
  ! the last dimension if it's MAF is used as an index into the arrays.

  subroutine JoinL2AUXQuantities ( key, name, quantity, l2auxDatabase, &
    & chunkNo, firstInstance, lastInstance, unpackOutput)

    ! Dummy arguments
    integer, intent(in) :: KEY     ! spec_args to decorate with the L2AUX index
    integer, intent(in) :: NAME    ! of the l2aux command
    type (VectorValue_T), intent(in) :: quantity
!   integer, intent(in) :: quantityNo
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    integer, intent(in), optional :: firstInstance, lastInstance
    logical, intent(in), optional :: unpackOutput
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2aux data.

    ! Local variables

    integer ::                              status, profNo
    integer ::                              useFirstInstance, useLastInstance, &
    &                                          noOutputInstances
    type (L2AUXData_T) ::                   newL2AUX
    type (L2AUXData_T), pointer ::          thisL2AUX
    logical ::                              l2auxDataIsNew, useUnpackOutput
    integer, dimension(3) ::                dimensionFamilies, dimensionSizes
    integer ::                              auxFamily     ! Channel or Frequency
    integer ::                              dimensionIndex,channel,surf,prof, &
    &                                         noMAFs,index, InstanceNo
    integer ::                              firstProfile,lastProfile
    real(r8), dimension(:,:), pointer ::    values
    INTEGER                              :: IERR

    ! Executable code

    if ( toggle(gen) ) call trace_begin ( "JoinL2AUXQuantities", key )

    ! If this is the first chunk, we have to setup the l2aux quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.

    l2auxDataIsNew = (.NOT. associated(l2auxDatabase))
    if ( .NOT. l2auxDataIsNew ) then
      index = decoration(key)
      l2auxDataIsNew = (index>=0)
    end if

    ! Work out what to do with the first and last Instance information
    
    if ( PRESENT(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( PRESENT(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    noOutputInstances = useLastInstance-useFirstInstance+1
    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    ! Sort out the unpackOutput flag

    if ( PRESENT(unpackOutput) ) then
      useUnpackOutput = unpackOutput
    else
      useUnpackOutput = .FALSE.
    end if

    if ( useUnpackOutput.AND.(.NOT.quantity%template%minorFrame) ) call MLSMessage ( &
      & MLSMSG_Error ,ModuleName, "Can only use the unpack output flag"// &
      & " for minor frame quantities")

    ! Now if this is a new l2aux quantity, we need to setup an l2aux data type
    ! for it.

    if ( l2auxDataIsNew ) then

      ! If the quantity is a minor frame quantity, then we deal with it 
      ! as such.  Otherwise we output it as a geodAngle based quantity

      if ( (quantity%template%noChans/=1) .AND. &
        & (quantity%template%frequencyCoordinate == L_None) ) &
        & CALL MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Quantity has multiple channels but no frequency coordinate" )

      select case (quantity%template%frequencyCoordinate)
      case ( L_InstrumentChannel )
        auxFamily = L2AUXDim_Channel
      case ( L_IntermediateFrequency )
        auxFamily = L2AUXDim_IntermediateFrequency
      case ( L_USBFrequency )
        auxFamily = L2AUXDim_USBFrequency
      case ( L_LSBFrequency )
        auxFamily = L2AUXDim_LSBFrequency
      case ( L_None ) ! OK to do nothing
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Unrecognised frequency coordinate" )
      end select

      if ( quantity%template%minorFrame ) then
        ! For minor frame quantities, the dimensions are:
        ! ([frequency or channel],MIF,MAF)
        !
        if ( useUnpackOutput ) then
          noMAFs = quantity%template%mafIndex(quantity%template%noInstances)
        else
          noMAFs = quantity%template%noInstances
        end if

        if ( quantity%template%frequencyCoordinate==L_None ) then
          dimensionFamilies = (/L2AUXDim_None, L2AUXDim_MIF, L2AUXDim_MAF/)
          dimensionSizes = (/1, quantity%template%noSurfs, noMAFs/)
        else
          dimensionFamilies = (/auxFamily, L2AUXDim_MIF, L2AUXDim_MAF/)
          dimensionSizes = (/quantity%template%noChans, quantity%template%noSurfs, &
            & noMAFs/)
        end if
      else
        ! Not a minor frame quantity; for non minor frame l2aux quantities
        ! our ability to output them will probably increase, but at the
        ! moment, I can't really forsee what form they may take.

        ! For the moment (ie. v0.1) I'm going to be restrictive and only
        ! allow quantities with no vertical coordinate.  This may and
        ! probably will change in later versions, leading to more L2AUXDim
        ! paramters etc., to parallel those from type t_verticalCoordinate
        ! in Init_Tables_Module.

        if ( quantity%template%verticalCoordinate /= l_None ) &
          & CALL MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Cannot currently output L2AUX quantities with obscure "// &
          & "vertical coordinates, sorry!" )

        if ( quantity%template%frequencyCoordinate==L_None ) then
          dimensionFamilies = (/L2AUXDim_geodAngle, L2AUXDim_None, &
            & L2AUXDim_None/)
          dimensionSizes = (/quantity%template%noInstances, 1, 1/)
        else
          dimensionFamilies = (/auxFamily, L2AUXDim_geodAngle, &
            & L2AUXDim_None/)
          dimensionSizes = (/quantity%template%noChans, quantity%template%noInstances, 1/)
        end if
      end if

      ! Now we setup the new quantity

      call SetupNewL2AUXRecord ( dimensionFamilies, dimensionSizes, newL2AUX )

      ! Setup the standard `vertical' and `channel' dimensions

      do dimensionIndex = 1, newL2aux%noDimensionsUsed-1
        select case ( dimensionFamilies(dimensionIndex) )
        case ( L2AUXDim_None )
          newL2AUX%dimensions(dimensionIndex)%values = 1
        case ( L2AUXDim_Channel )
          do channel = 1,quantity%template%noChans
            newL2AUX%dimensions(dimensionIndex)%values(channel) = channel
          end do
        case ( L2AUXDim_IntermediateFrequency )
          newL2AUX%dimensions(dimensionIndex)%values = quantity%template%frequencies
        case ( L2AUXDim_USBFrequency )
          newL2AUX%dimensions(dimensionIndex)%values = quantity%template%frequencies
        case ( L2AUXDim_LSBFrequency )
          newL2AUX%dimensions(dimensionIndex)%values = quantity%template%frequencies
        case ( L2AUXDim_MIF )
          do surf = 1, quantity%template%noSurfs
            newL2AUX%dimensions(dimensionIndex)%values(surf) = surf
          end do
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Obscure error with coordinates in l2aux data" )
        end select
        ! The error message here is rather vaugue.  The issue is that
        ! both MAF and geodAngle should only occur for the `last' dimension
        ! which our loop is explicity avoiding.
      end do ! The `last' dimension is dealt with later on.

      ! Add this l2aux to the database
      index = AddL2AUXToDatabase ( l2auxDatabase, newL2AUX )

      ! Setup the pointer and the index to be used later
      call decorate ( key, -index ) ! Remember where it is
      thisL2AUX => l2auxDatabase(index)

      ! Finally destroy the information in the interim, non database l2aux
!     call DestroyL2AUXContents ( newL2AUX )
      ! No! Don't do this! It deallocates all the pointers we just
      ! copied into the database!
    else
      ! Setup the index and pointer
      thisL2AUX => l2auxDatabase(-index)

      ! Expand this l2aux along the `last' dimension to take up the new
      ! information.

      ! ??? The noMAFs computation isn't consistent with the initial case ???
      if ( useUnpackOutput ) then
        noMAFs = thisL2AUX%dimensions(thisL2AUX%noDimensionsUsed)%noValues+&
          & quantity%template%noInstances
      else
        noMAFs = quantity%template%mafCounter(quantity%template%noInstances)
      end if

      ! ??? WVS would like to make this work as in the L2GP case:  The
      ! "setup..." routine allocates zero size in the MAF's direction,
      ! Then "expand..." is called in either case (but does no copying when
      ! called immediately after the initial one).
      call ExpandL2AUXDataInPlace ( thisL2AUX, noMAFs )
    end if

    ! Now we are ready to fill up the l2aux quantity with the new data.
    ! We do this differently depending on whether we want packed or unpacked
    ! data.

    thisL2AUX%name = name

    if ( .NOT. useUnpackOutput ) then
!     lastProfile = thisL2AUX%dimensions(thisL2AUX%noDimensionsUsed)%noValues
!     firstProfile = lastProfile-noOutputInstances+1
    end if

!    vectorQuantity = GetVectorQuantity ( vector, quantity%template%name, &
!      & quantityIsName=.TRUE. )

!    thisL2AUX%values = &
!      & vectorQuantity%values(useFirstInstance:useLastInstance,:,:)
    CALL unsqueeze(quantity%values(:, useFirstInstance:useLastInstance), &
    & thisL2AUX%values, IERR)

    if ( toggle(gen) ) call trace_end ( "JoinL2AUXQuantities" )

  end subroutine JoinL2AUXQuantities

  ! ---------------------------------------------------------------------------

!=============================== unsqueeze ==========================
SUBROUTINE unsqueeze(source, sink, IERR)
!=============================== unsqueeze ==========================
    ! takes a rank2 object source and returns a rank3 object sink
    ! source(1..n1*n2, 1..n3) -> sink(1..n1, 1..n2, 1..n3)
    ! unless it can't--then it returns iERR /= 0
    ! One reason it may fail: shape of sink too small
    !
    ! Assuming that shape(sink) = {n1, n2, n3}
    !     =>      shape(source) = {m1, m2}
    ! then we must further assume (else set IERR)
    ! m1 <= n1*n2
    ! m2 <= n3
    !
    ! (A future improvement might take as optional arguments
    !  integer arrays source_shape, sink_shape, 
    !  or else shape-params n1, n2, m1)
    !--------Argument--------!
    REAL(r8), DIMENSION(:,:), INTENT(IN)     :: source
    REAL(r8), DIMENSION(:,:,:), INTENT(OUT)  :: sink
    INTEGER, INTENT(OUT)                     :: IERR
 
    !----------Local vars----------!
    ! Error codes
    INTEGER               :: m1_is_zero = 1
    INTEGER               :: m2_is_zero = 2
    INTEGER               :: m1_too_big = 4
    INTEGER               :: m2_too_big = 5
    
    INTEGER, DIMENSION(4) :: source_shape
    INTEGER, DIMENSION(4) :: sink_shape
    INTEGER::i,icode,offset
    !----------Executable part----------!
    source_shape(1:2) = shape(source)
    sink_shape(1:3) = shape(sink)

    IF (source_shape(1) == 0) THEN
    	IERR = m1_is_zero
        RETURN
    ELSEIF (source_shape(2) == 0) THEN
    	IERR = m2_is_zero
        RETURN
    ELSEIF (sink_shape(1)*sink_shape(2) < source_shape(1)) THEN
    	IERR = m1_too_big
        RETURN
    ELSEIF (sink_shape(3) < source_shape(2)) THEN
    	IERR = m2_too_big
        RETURN
    ELSE
    	IERR = 0
    ENDIF

    sink = reshape(source, sink_shape(1:3))
    
END SUBROUTINE unsqueeze

!=============================================================================
end module Join
!=============================================================================

!
! $Log$
! Revision 2.15  2001/03/05 01:19:45  livesey
! Removed a print statement
!
! Revision 2.14  2001/03/05 01:01:12  livesey
! Bug fix, now uses GetVectorQtyFromTemplateIndex
!
! Revision 2.13  2001/03/01 18:38:27  livesey
! Fixed bug with verticalCoordinate==l_Zeta
!
! Revision 2.12  2001/02/27 17:38:21  livesey
! Tidied things up, removed unnecessary arguments
!
! Revision 2.11  2001/02/27 00:50:31  livesey
! Added ability to Join verticalCoordinate=l_zeta quantities into l2gp entities.
!
! Revision 2.10  2001/02/16 00:50:17  livesey
! Added error to avoid confusion with L2GP in ReadApriori section
!
! Revision 2.9  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.8  2001/02/09 18:01:46  livesey
! Various further updates, set default values for status and quality
!
! Revision 2.7  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.6  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.5  2000/11/16 02:19:01  vsnyder
! Implement timing.
!
! Revision 2.4  2000/11/13 23:02:21  pwagner
! Adapted for rank2 vectorsModule
!
! Revision 2.3  2000/10/05 16:37:19  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:34:35  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!




