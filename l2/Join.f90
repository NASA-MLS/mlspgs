! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Join                     ! Join together chunk based data.
!=============================================================================

  ! This module performs the 'join' task in the MLS level 2 software.

  use INIT_TABLES_MODULE, only: F_COMPAREOVERLAPS, F_FILE, F_OUTPUTOVERLAPS, &
    & F_PRECISION, F_PREFIXSIGNAL, F_SOURCE, F_SDNAME, F_SWATH, FIELD_FIRST, &
    & FIELD_LAST
  use INIT_TABLES_MODULE, only: L_PRESSURE, L_NONE, &
    & L_TRUE, L_ZETA, S_L2AUX, S_L2GP, S_TIME
  use Intrinsic, ONLY: FIELD_INDICES, L_NONE, L_CHANNEL, L_GEODANGLE, &
    & L_INTERMEDIATEFREQUENCY, L_LSBFREQUENCY, L_MAF, L_MIF, L_USBFREQUENCY
  use L2AUXData, only: AddL2AUXToDatabase, ExpandL2AUXDataInPlace, &
    & L2AUXData_T, L2AUXRank, SetupNewL2AUXRecord
  use L2GPData, only: AddL2GPToDatabase, ExpandL2GPDataInPlace, &
    & L2GPData_T, SetupNewL2GPRecord
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: MLSChunk_T, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
    & MLSMSG_Allocate, MLSMSG_Deallocate
  use MLSSignals_M, only: GETSIGNALNAME
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: DISPLAY_STRING, GET_STRING
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_STRING
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
    & SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_SET_ONE
  use VectorsModule, only: GetVectorQuantity, GetVectorQtyByTemplateIndex, &
    & ValidateVectorQuantity, Vector_T, VectorValue_T, DUMP

  implicit none
  private
  public :: MLSL2Join

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
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

  subroutine MLSL2Join ( root, vectors, l2gpDatabase, l2auxDatabase, &
    & chunkNo, chunks )

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the JOIN section in the AST
    type (Vector_T), dimension(:), pointer :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks

    ! Local variables
    logical :: COMPAREOVERLAPS
    integer :: FIELD                    ! Subtree index of "field" node
    integer :: FIELD_INDEX              ! F_..., see Init_Tables_Module
    integer :: HDFNAMEINDEX             ! Name of swath/sd
    integer :: GSON                     ! Son of Key
    integer :: KEY                      ! Index of an L2GP or L2AUX tree
    integer :: KEYNO                    ! Index of subtree of KEY
    integer :: MLSCFLine
    logical :: OutputOverlaps
    integer :: NAME                     ! Sub-rosa index of name of L2GP or L2AUX
    integer :: SON                      ! A son of ROOT
    integer :: SOURCE                   ! Index in AST
    integer :: STATUS                   ! Flag
    integer :: VALUE                    ! Value of a field
    integer :: VECTORINDEX              ! Index for vector to join
    integer :: QUANTITYINDEX            ! ind in qty tmpl database, not vector
    logical :: PREFIXSIGNAL             ! Prefix (i.e. make) the sd name the signal
    integer :: PRECVECTORINDEX          ! Index for precision vector
    integer :: PRECQTYINDEX             ! Index for precision qty (in database not vector)
    logical :: TIMING

    real :: T1, T2     ! for timing

    character(len=132) :: HDFNAME          ! Name for swath/sd
    logical :: GOT_FIELD(field_first:field_last)
    type (VectorValue_T), pointer :: Quantity
    type (VectorValue_T), pointer :: PrecisionQuantity

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
      select case( get_spec_id(key) )
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
      compareOverlaps = .FALSE.
      outputOverlaps = .FALSE.
      hdfNameIndex=name
      prefixSignal = .false.

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
        case ( f_precision )
          source = subtree(2,gson) ! required to be an n_dot vertex
          precVectorIndex = decoration(decoration(subtree(1,source)))
          precQtyIndex = decoration(decoration(decoration(subtree(2,source))))
        case ( f_prefixSignal )
          prefixSignal= value == l_true
        case ( f_compareoverlaps )
          compareOverlaps = value == l_true
        case ( f_outputoverlaps )
          outputOverlaps = value == l_true 
        case (f_swath)
          hdfNameIndex = sub_rosa(subtree(2,gson))
        case (f_sdName)
          hdfNameIndex = sub_rosa(subtree(2,gson))
       case ( f_file)
          call announce_error(key,NotAllowed,field_index)
        case default ! Can't get here if tree_checker worked properly
        end select
      end do

      if ( error > 0 ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Errors in configuration prevent proceeding" )

      select case ( get_spec_id(key) )
      case ( s_l2gp, s_l2aux ) ! ------------- L2GP and L2AUX Data     ------------
        quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex),quantityIndex)
        if ( got_field ( f_precision ) ) then
          precisionQuantity => &
            & GetVectorQtyByTemplateIndex(vectors(precVectorIndex),precQtyIndex)
          if ( quantity%template%id /= precisionQuantity%template%id ) &
            & call MLSMessage(MLSMSG_Error,ModuleName, &
            & 'Quantity and precision quantity do not match')
        else
          precisionQuantity => NULL()
        endif
        
        hdfName = ''
        if ( prefixSignal ) &
          & call GetSignalName ( quantity%template%signal, hdfName, &
          &   sideband=quantity%template%sideband )
        call Get_String( hdfNameIndex, hdfName(len_trim(hdfName)+1:), strip=.true. )
        hdfNameIndex = enter_terminal ( trim(hdfName), t_string )

        ! Now, depending on the properties of the source we deal with the
        ! vector quantity appropriately.
        if (ValidateVectorQuantity(quantity,coherent=.TRUE.,stacked=.TRUE.,regular=.TRUE.,&
          & verticalCoordinate=(/L_Pressure,L_Zeta,L_None/))) then
          ! Coherent, stacked, regular quantities on pressure surfaces, or
          ! with no vertical coordinate system go in l2gp files.
          if ( get_spec_id(key) /= s_l2gp ) call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2gp')
          call JoinL2GPQuantities ( key, hdfNameIndex, quantity, precisionQuantity, &
            & l2gpDatabase, chunkNo )
        else
          ! All others go in l2aux files.
          if ( get_spec_id(key) /= s_l2aux ) call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2aux')
          call JoinL2AUXQuantities ( key, hdfNameIndex, quantity, &
            & l2auxDatabase, chunkNo, chunks )
        endif

      case default ! Timing
      end select

    end do

    if ( toggle(gen) ) call trace_end ( "MLSL2Join" )
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

  subroutine JoinL2GPQuantities ( key, name, quantity, precision, l2gpDatabase, &
    & chunkNo, firstInstance, lastInstance )

    ! Dummy arguments
    integer, intent(in) :: KEY          ! spec_args to Decorate with the L2GP index
    integer, intent(in) :: NAME         ! For the swath
    type (VectorValue_T), intent(in) :: QUANTITY ! Vector quantity
    type (VectorValue_T), pointer :: PRECISION ! Optional vector quantity
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    integer, intent(in) :: CHUNKNO
    integer, intent(in), optional :: FIRSTINSTANCE, LASTINSTANCE
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2gp data.

    ! Local variables

    integer :: FreqNo, ProfNo, Status
    type (L2GPData_T) :: NewL2GP
    type (L2GPData_T), pointer :: ThisL2GP
    integer :: Index
    integer :: FirstProfile, LastProfile ! Profile range in the l2gp to output to
    integer :: NoSurfsInL2GP, NoFreqsInL2GP
    integer :: UseFirstInstance, UseLastInstance, NoOutputInstances
    logical :: L2gpDataIsNew
!   real(r8), dimension(:,:), pointer :: Values !??? Not used ???
    
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
    call Get_String( name, thisL2GP%name, strip=.true.)
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

    thisL2GP%l2gpValue(:,:,firstProfile:lastProfile) = &
         RESHAPE(quantity%values(:,useFirstInstance:useLastInstance),&
         (/MAX(thisL2GP%nFreqs,1),MAX(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    if (associated(precision)) &
      thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = &
      RESHAPE(precision%values(:,useFirstInstance:useLastInstance),&
         (/MAX(thisL2GP%nFreqs,1),MAX(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0 ! Later put something here
    thisL2GP%status(firstProfile:lastProfile)='G'
    thisL2GP%quality(firstProfile:lastProfile)=0.0

    if ( toggle(gen) ) call trace_end ( "JoinL2GPQuantities" )
  end subroutine JoinL2GPQuantities

  ! ----------------------------------------  JoinL2AUXQuantities  -----

  ! This subroutine is like the one above, except that the quantities it joins
  ! are destined to go in L2AUX quantities.

  subroutine JoinL2AUXQuantities ( key, name, quantity, l2auxDatabase, &
    & chunkNo, chunks, firstInstance, lastInstance )

    ! Dummy arguments
    integer, intent(in) :: KEY     ! spec_args to decorate with the L2AUX index
    integer, intent(in) :: NAME    ! for the sd
    type (VectorValue_T), intent(in) :: quantity
!   integer, intent(in) :: quantityNo
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    integer, intent(in), optional :: firstInstance, lastInstance
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2aux data.

    ! Local variables
    integer :: FIRSTMAF                 ! Index
    integer :: LASTMAF                  ! Index
    integer ::                           Status, ProfNo
    integer ::                           UseFirstInstance, UseLastInstance, &
    &                                    NoOutputInstances
    type (L2AUXData_T) ::                NewL2AUX
    type (L2AUXData_T), pointer ::       ThisL2AUX
    logical ::                           L2auxDataIsNew
    integer, dimension(3) ::             DimensionFamilies, DimensionSizes, DimensionStarts
    integer ::                           AuxFamily     ! Channel or Frequency
    integer ::                           DimensionIndex, Channel, Surf, Prof, &
    &                                    NoMAFs,index, InstanceNo
    integer ::                           FirstProfile, LastProfile
!   real(r8), dimension(:,:), pointer :: values !??? Not used ???
    integer ::                           IERR

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

    ! Now if this is a new l2aux quantity, we need to setup an l2aux data type
    ! for it.

    if ( l2auxDataIsNew ) then

      ! If the quantity is a minor frame quantity, then we deal with it 
      ! as such.  Otherwise we output it as a geodAngle based quantity

      if ( (quantity%template%noChans/=1) .AND. &
        & (quantity%template%frequencyCoordinate == L_None) ) &
        & CALL MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Quantity has multiple channels but no frequency coordinate" )

      auxFamily=quantity%template%frequencyCoordinate

      if ( quantity%template%minorFrame ) then
        ! For minor frame quantities, the dimensions are:
        ! ([frequency or channel],MIF,MAF)
        !
        ! For minor frame quantities, we're going to allocate them to the full
        ! size from the start, rather than expand them, as this will be more
        ! efficient
        firstMAF = minval ( chunks%firstMAFIndex )
        lastMAF = maxval ( chunks%lastMAFIndex )
        noMAFs = lastMAF - firstMAF + 1
        ! THINK HERE ABOT RUNS THAT DON'T START AT THE BEGINNING !???????
        ! MAY NEED TO CHANGE ALLOCATE
        if ( quantity%template%frequencyCoordinate==L_None ) then
          dimensionFamilies = (/L_None, L_MIF, L_MAF/)
          dimensionSizes = (/1, quantity%template%noSurfs, noMAFs/)
          dimensionStarts = (/1, quantity%template%noSurfs, firstMAF /)
        else
          dimensionFamilies = (/auxFamily, L_MIF, L_MAF/)
          dimensionSizes = (/quantity%template%noChans, quantity%template%noSurfs, &
            & noMAFs/)
          dimensionStarts = (/1, 1, firstMAF /)
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
          dimensionFamilies = (/L_geodAngle, L_None, &
            & L_None/)
          dimensionSizes = (/quantity%template%noInstances, 1, 1/)
          dimensionStarts = (/1, 1, 1/)
        else
          dimensionFamilies = (/auxFamily, L_geodAngle, &
            & L_None/)
          dimensionSizes = (/quantity%template%noChans, quantity%template%noInstances, 1/)
          dimensionStarts = (/1, 1, 1/)
        end if
      end if

      ! Now we setup the new quantity
      call SetupNewL2AUXRecord ( dimensionFamilies, dimensionSizes, &
        & dimensionStarts, newL2AUX )
      newL2AUX%minorFrame=quantity%template%minorFrame
      newL2AUX%instrumentModule=quantity%template%instrumentModule

      ! Setup the standard `vertical' and `channel' dimensions

      do dimensionIndex = 1, L2AUXRank
        select case ( dimensionFamilies(dimensionIndex) )
        case ( L_None )          ! Do nothing
        case ( L_Channel )
          do channel = 1,quantity%template%noChans
            newL2AUX%dimensions(dimensionIndex)%values(channel) = channel
          end do
        case ( L_IntermediateFrequency, l_USBFrequency, L_LSBFrequency )
          newL2AUX%dimensions(dimensionIndex)%values = quantity%template%frequencies
        case ( L_MIF )
          do surf = 1, quantity%template%noSurfs
            newL2AUX%dimensions(dimensionIndex)%values(surf) = surf
          end do
        case default                    ! Ignore horizontal dimensions
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

    else
      ! Setup the index and pointer
      thisL2AUX => l2auxDatabase(-index)

      ! Expand this l2aux along the `last' dimension to take up the new
      ! information.

      ! ??? The noMAFs computation isn't consistent with the initial case ???
      noMAFs = quantity%template%mafCounter( &
        & quantity%template%noInstances-quantity%template%noInstancesUpperOverlap)

      ! ??? WVS would like to make this work as in the L2GP case:  The
      ! "setup..." routine allocates zero size in the MAF's direction,
      ! Then "expand..." is called in either case (but does no copying when
      ! called immediately after the initial one).
      ! For minor frame quantities we don't need to expand, as they're created
      ! at full size from the start.
      if (.not. quantity%template%minorFrame) &
        & call ExpandL2AUXDataInPlace ( thisL2AUX, noMAFs )
    end if

    ! Now we are ready to fill up the l2aux quantity with the new data.
    thisL2AUX%name = name

    if (quantity%template%minorFrame) then
      lastProfile = quantity%template%mafIndex(quantity%template%noInstances)
    else
      lastProfile = thisL2AUX%dimensions(L2AUXRank)%noValues
    endif
    firstProfile = lastProfile-noOutputInstances+1
    
    select case (thisL2AUX%dimensions(L2AUXRank)%dimensionFamily)
    case ( L_GeodAngle )
      thisL2AUX%dimensions(L2AUXRank)%values(firstProfile:lastProfile)=&
        & quantity%template%phi(1,useFirstInstance:useLastInstance)
    case ( L_MAF )
      thisL2AUX%dimensions(L2AUXRank)%values(firstProfile:lastProfile)=&
        & quantity%template%mafCounter(useFirstInstance:useLastInstance)
    case default
    end select
    
    thisL2AUX%values(:,:,firstProfile:lastProfile) = &
      & reshape(quantity%values(:,useFirstInstance:useLastInstance), &
      &   (/ thisL2AUX%dimensions(1)%noValues, &
      &      thisL2AUX%dimensions(2)%noValues, &
      &      lastProfile-firstProfile+1/) )
    
    if ( toggle(gen) ) call trace_end ( "JoinL2AUXQuantities" )

  end subroutine JoinL2AUXQuantities

end module Join

!
! $Log$
! Revision 2.37  2001/05/12 00:18:17  livesey
! Tidied up array bounds for L2AUX/MAF.
!
! Revision 2.36  2001/05/10 16:31:24  livesey
! Added prefix signal option for swath/sd name
!
! Revision 2.35  2001/05/08 23:25:32  livesey
! Added the precision stuff for l2gp's
!
! Revision 2.34  2001/05/08 21:51:02  livesey
! Removed some old xStar, yStar, kStar stuff.
!
! Revision 2.33  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.32  2001/05/02 22:22:43  pwagner
! Removed SDPToolkit use
!
! Revision 2.31  2001/04/28 01:30:14  livesey
! Basically gone back to an earlier version.  As l2pc's now output
! directly as matrices there is no need for Join to think about them.
!
! Revision 2.30  2001/04/27 21:52:39  livesey
! Removed l2pc stuff
!
! Revision 2.29  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.28  2001/04/26 15:59:30  livesey
! Tidied up uses
!
! Revision 2.27  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.26  2001/04/26 00:07:16  livesey
! Insulate vector is gone
!
! Revision 2.25  2001/04/25 21:54:22  livesey
! Added candol2pc flag
!
! Revision 2.24  2001/04/25 21:51:46  livesey
! Tidied up Join for l2pcs
!
! Revision 2.23  2001/04/25 20:33:28  livesey
! Minor improvements to Join l2pc stuff
!
! Revision 2.22  2001/04/24 20:20:27  livesey
! L2PC moved to lib and word bin dropped from types etc.
!
! Revision 2.21  2001/04/24 20:04:54  livesey
! Added l2pc joining
!
! Revision 2.20  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.19  2001/03/15 23:26:56  livesey
! Avoid calling ExpandL2AUXQuantitiesInPlace for minor frame quantities.
! Really saves on memory thrashing.
!
! Revision 2.18  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.17  2001/03/06 22:40:41  livesey
! New L2AUX stuff
!
! Revision 2.16  2001/03/05 20:46:41  livesey
! Removed a debugging statement left behind
!
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

