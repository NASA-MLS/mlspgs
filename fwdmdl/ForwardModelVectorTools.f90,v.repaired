! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ForwardModelVectorTools          ! Tools for vectors in forward models

  ! This module contains routines needed to help a forward model get
  ! hold of the quantities it needs.

  implicit NONE

  private

  public :: GetQuantityForForwardModel, GetQtyStuffForForwardModel

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------  GetQuantityForForwardModel  -----
  function GetQuantityForForwardModel ( vector, otherVector, quantityType, &
    & molecule, instrumentModule, supportedInstrumentModule, radiometer, &
    & reflector, signal, sideband, &
    & molIndex, config, foundInFirst, wasSpecific, noError, matchQty, Frq )

    ! This function is in many senses like GetVectorQuantityByType, (to
    ! which it can revert), except that given a forwardModelConfig_T in
    ! config, and possibly an index into the molecules array, it can
    ! use the specificQuantities stuff in config to identify exactly
    ! the right quantity.
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T
    use INTRINSIC, only: LIT_INDICES
    use INTRINSIC, only: L_VMR
    use MANIPULATEVECTORQUANTITIES, only: DOHGRIDSMATCH, DOVGRIDSMATCH
    use MLSKINDS, only: R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSFINDS, only: FINDFIRST
    use MLSSIGNALS_M, only: GETRADIOMETERNAME, GETSIGNALNAME, SIGNALS
    use MOLECULES, only: ISEXTINCTION
    use QUANTITYTEMPLATES, only: QUANTITYTEMPLATE_T
    use STRING_TABLE, only: GET_STRING
    use VECTORSMODULE, only: GETVECTORQUANTITYBYTYPE, VECTOR_T, VECTORVALUE_T

    ! Dummy arguments
    type (Vector_T), target :: VECTOR ! First vector to look in
    type (Vector_T), target, optional :: OTHERVECTOR ! Second vector to look in
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in),  optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in),  optional :: INSTRUMENTMODULE ! Instrument module index
    integer, intent(in),  optional :: SUPPORTEDINSTRUMENTMODULE ! Another instrument module index
    integer, intent(in),  optional :: RADIOMETER   ! Radiometer index
    integer, intent(in),  optional :: REFLECTOR    ! Reflector literal
    integer, intent(in),  optional :: SIGNAL       ! Signal index
    integer, intent(in),  optional :: SIDEBAND     ! -1, 0, +1
    type (ForwardModelConfig_T), intent(in), optional :: CONFIG ! fwmConfig
    integer, intent(in),  optional :: MOLINDEX     ! Index into the molecules array
    logical, intent(out), optional :: FOUNDINFIRST ! Set if found in first vector
    logical, intent(out), optional :: WASSPECIFIC  ! Set if listed as specific quantity
    logical, intent(in),  optional :: NOERROR      ! Don't give error if not found
    type (VectorValue_T), intent(in), optional :: MATCHQTY ! Result must match this
    real(r8), intent(in), optional :: Frq          ! Frequency
    ! Result
    type (VectorValue_T), pointer :: GetQuantityForForwardModel

    ! Local type
    type Stuff_t ! a vector_t pointer and a real(r8) pointer
      type(vector_t), pointer :: V => null() ! Vector(1) or OtherVector(2)
      integer, pointer :: Match(:) => null() ! Flags for V
    end type

    ! Local variables
    logical :: MyNoError                ! Copy of no error
    real(r8) :: FrqDiff

    ! -1 for no match.  Channel number of best match if Frq is present, else zero.
    integer, dimension(:), pointer :: MATCH   ! One of stuff(:)%match

    real(r8) :: MyFrq                   ! Copy of Frq, or a fiction if not present

    integer :: BestMatch                ! Channel index of best match
    integer :: MolEntry                 ! How many similar molecules in list?
    integer :: NoFound                  ! Number of matches found so far.
    integer :: NoVectors                ! Number of vectors we've been given
    integer :: Quantity                 ! Loop counter
    integer :: ThisMolecule             ! A molecule
    integer :: VectorIndex              ! Loop counter
    type (Vector_T), pointer :: V       ! A vector
    type (QuantityTemplate_T), pointer :: QT ! A quantity template
    character(len=127) :: MSG
    logical :: UseGetQuantityByType

    ! pointers to vectors and Match arrays
    type(stuff_t) :: stuff(2)

    ! Executable code

    ! Do a sanity check, get default values
    if ( present ( molIndex ) .and. .not. present ( config ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot have molIndex in GetQuantityForForwardModel without config" )
    myNoError = .false.
    if ( present ( noError ) ) myNoError = noError
    if ( present ( wasSpecific ) ) wasSpecific = .false.
    if ( present ( foundInFirst ) ) foundInFirst = .false.

    myFrq = huge(0.0_r8)
    if ( present(frq) ) myFrq = frq
    nullify ( GetQuantityForForwardModel )

    ! First see if we can simply revert to the simpler GetVectorQuantityByType
    useGetQuantityByType = .true.
    if ( present ( config ) ) then
      if ( associated ( config%specificQuantities ) .or. present(molIndex) ) &
        & useGetQuantityByType = .false.
    end if

    ! If we can revert to the simpler GetVectorQuantityByType then do so.
    if ( useGetQuantityByType ) then
      GetQuantityForForwardModel => GetVectorQuantityByType ( vector, otherVector, &
        & quantityType, molecule, instrumentModule, supportedInstrumentModule, &
        & radiometer, reflector, signal, sideband, &
        & foundInFirst, noError )
      return
    end if

    ! OK, looks like we have to do it the more complicated way.
    stuff(1)%v => vector
    if ( present(otherVector) ) stuff(2)%v => otherVector
    noVectors = merge(2,1,present(otherVector))

    do vectorIndex = 1, noVectors
      call allocate_test ( stuff(vectorIndex)%match, &
        & stuff(vectorIndex)%v%template%noQuantities, &
        & 'stuff%match(vectorIndex)', ModuleName, fill=-1 )
    end do

    ! Now loop over the one or two vectors and work out which quantities might match.
    do vectorIndex = 1, noVectors
      v => stuff(vectorIndex)%v
      match => stuff(vectorIndex)%match
      do quantity = 1, v%template%noQuantities
        thisMolecule = 0
        qt => v%quantities(quantity)%template
        ! Now go through the quantities and see if they match
        if ( quantityType /= qt%quantityType ) cycle
        if ( present(molecule) ) then
          if ( molecule /= qt%molecule ) cycle
          thisMolecule = molecule
        end if
        if ( present(molIndex) ) then
          if ( config%molecules ( molIndex ) /= qt%molecule ) cycle
          thisMolecule = config%molecules ( molIndex )
        end if
        if ( present(instrumentModule ) ) then
          if ( instrumentModule /= qt%instrumentModule ) cycle
        end if
        if ( present(radiometer) ) then
          ! We can be a little lenient here in the case of vmrs
          if ( quantityType == l_vmr ) then
            if ( radiometer /= qt%radiometer .and. &
               & isExtinction(thisMolecule) ) cycle
          else
            if ( radiometer /= qt%radiometer ) cycle
          end if
        end if
        if ( present(reflector) ) then
          if ( reflector /= qt%reflector ) cycle
        end if
        if ( present(sideband) ) then
          if ( sideband /= qt%sideband ) cycle
        end if
        if ( present(signal) ) then
          if ( signal /= qt%signal ) cycle
        end if
        if ( present(matchQty) ) then
          if ( .not. DoVGridsMatch ( v%quantities(quantity), matchQty ) ) cycle
          if ( .not. DoHGridsMatch ( v%quantities(quantity), matchQty, &
            & spacingOnly=.true. ) ) cycle
        end if
        if ( qt%signal <= 0 .or. .not. present(frq) ) then
          match ( quantity ) = 0
        else
          frqDiff = huge(0.0_r8)
          do bestMatch = 1, size(signals(qt%signal)%frequencies)
            if ( signals(qt%signal)%channels(bestMatch) ) then
              if ( abs(signals(qt%signal)%frequencies(bestMatch)-myFrq) < frqDiff ) then
                frqDiff = abs(signals(qt%signal)%frequencies(bestMatch)-myFrq)
                if ( frqDiff <= signals(qt%signal)%widths(match(quantity)) ) &
                  & match ( quantity ) = bestMatch
              end if
            end if
          end do
        end if
      end do ! Quantity                 ! End loop over the quantities
    end do   ! VectorIndex              ! End loop over the one or two vectors

    ! Now decide exactly which ones we want
    if ( present ( molIndex ) ) then
      ! Some special thought for this case.
      ! Work out how many of the same molecules preceed or include this one.
      molEntry = count ( config%molecules ( 1:molIndex ) == &
        & config%molecules ( molIndex ) )
    else
      molEntry = 0
    end if

    noFound = 0

    ! If Frq is present, find the quantity for a signal whose center frequency
    ! is nearest to the desired one, and verify that the desired frequency is
    ! within the same channel.
    if ( present(frq) ) then

    ! First to see if any of the matches are on our 'specificQuantity'
    ! list, if so match them.  Though pay special attention in the molIndex case.
    else if ( present(config) ) then
      if ( associated ( config%specificQuantities ) ) then
        specificVectorLoop: do vectorIndex = 1, noVectors
          v => stuff(vectorIndex)%v
          match => stuff(vectorIndex)%match
          do quantity = 1, v%template%noQuantities
            if ( match ( quantity ) >= 0 .and. &
              & any ( v%template%quantities(quantity) == config%specificQuantities ) ) then
              noFound = noFound + 1
              if ( molEntry == 0 .or. molEntry == noFound ) then
                GetQuantityForForwardModel => v%quantities(quantity)
                if ( present ( foundInFirst ) ) foundInFirst = ( vectorIndex == 1 )
                if ( present ( wasSpecific ) ) wasSpecific = .true.
                exit specificVectorLoop
              end if
            end if
          end do
        end do specificVectorLoop
      end if
    end if

    if ( molEntry > 1 .and. noFound < molEntry ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unresolved ambiguity in molecule list' )

    ! Otherwise just get it from one or other vector
    if ( noFound == 0 ) then
      nonSpecificVectorLoop: do vectorIndex = 1, noVectors
        v => stuff(vectorIndex)%v
        match => stuff(vectorIndex)%match
        quantity = FindFirst ( match >= 0 )
        if ( quantity /= 0 ) then
          GetQuantityForForwardModel => v%quantities(quantity)
          noFound = 1
          if ( present ( foundInFirst ) ) foundInFirst = ( vectorIndex == 1 )
          exit nonSpecificVectorLoop
        end if
      end do nonSpecificVectorLoop
    end if

    ! Now, if appropriate, print out an explanatory error message.
    if ( noFound == 0 .and. .not. myNoError ) then
      msg = 'There is no quantity in vector '
      if ( vector%name /= 0 ) then
        call get_string ( vector%name, msg(len_trim(msg)+2:) )
      else
        msg(len_trim(msg)+2:) = '[unnamed]'
      end if
      msg = trim(msg) // ' that has type'
      call get_string ( lit_indices(quantityType), msg(len_trim(msg)+2:) )

      if ( present ( molecule ) ) then
        msg = trim(msg) // ' for molecule'
        call get_string ( lit_indices(molecule), msg(len_trim(msg)+2:))
      end if
      if ( present ( molIndex ) ) then
        msg = trim(msg) // ' for molecule'
        call get_string ( lit_indices(config%molecules(molIndex)), msg(len_trim(msg)+2:))
      end if

      if ( present ( radiometer ) ) then
        msg = trim(msg) // ' for radiometer'
        call getRadiometerName ( radiometer, msg(len_trim(msg)+2:))
      end if

      if ( present ( reflector ) ) then
        msg = trim(msg) // ' for reflector'
        call get_string ( lit_indices(reflector), msg(len_trim(msg)+2:))
      end if

      if ( present ( instrumentModule ) ) then
        msg = trim(msg) // ' for instrument module'
        call get_string ( lit_indices(instrumentModule), msg(len_trim(msg)+2:))
      end if

      if ( present ( signal ) ) then
        msg = trim(msg) // ' for signal'
        call GetSignalName ( signal, msg(len_trim(msg)+2:), sideband=sideband )
      end if

      call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )
    end if

    call Deallocate_test ( stuff(1)%match, 'stuff%match(vectorIndex)', ModuleName )
    call Deallocate_test ( stuff(2)%match, 'stuff%match(vectorIndex)', ModuleName )

  end function GetQuantityForForwardModel


  ! ---------------------------------  GetQtyStuffForForwardModel  -----
  function GetQtyStuffForForwardModel ( vector, otherVector, quantityType, &
    & molecule, instrumentModule, supportedInstrumentModule, radiometer, reflector, signal, sideband, &
    & molIndex, config, noError, matchQty, Frq )

    ! This function is in many senses like GetVectorQuantityByType, (to
    ! which it can revert), except that given a forwardModelConfig_T in
    ! config, and possibly an index into the molecules array, it can
    ! use the specificQuantities stuff in config to identify exactly
    ! the right quantity.
    use FORWARDMODELCONFIG, only: FORWARDMODELCONFIG_T, QTYSTUFF_T
    use MLSKINDS, only: R8
    use VECTORSMODULE, only: VECTOR_T, VECTORVALUE_T

    ! Dummy arguments
    type (Vector_T), target :: VECTOR ! First vector to look in
    type (Vector_T), target, optional :: OTHERVECTOR ! Second vector to look in
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in),  optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in),  optional :: INSTRUMENTMODULE ! Instrument module index
    integer, intent(in),  optional :: SUPPORTEDINSTRUMENTMODULE ! Another instrument module index
    integer, intent(in),  optional :: RADIOMETER   ! Radiometer index
    integer, intent(in),  optional :: REFLECTOR    ! Reflector literal
    integer, intent(in),  optional :: SIGNAL       ! Signal index
    integer, intent(in),  optional :: SIDEBAND     ! -1, 0, +1
    type (ForwardModelConfig_T), intent(in), optional :: CONFIG ! fwmConfig
    integer, intent(in),  optional :: MOLINDEX     ! Index into the molecules array
    logical, intent(in),  optional :: NOERROR      ! Don't give error if not found
    type (VectorValue_T), intent(in), optional :: MATCHQTY ! Result must match this
    real(r8), intent(in), optional :: Frq          ! Frequency
    ! Result
    type (QtyStuff_t) :: GetQtyStuffForForwardModel

    GetQtyStuffForForwardModel%qty => GetQuantityForForwardModel ( vector, &
      & otherVector, quantityType, molecule, instrumentModule, supportedInstrumentModule, &
      & radiometer, reflector, signal, sideband, molIndex, config, &
      & getQtyStuffForForwardModel%foundInFirst, &
      & getQtyStuffForForwardModel%wasSpecific, noError, matchQty, Frq )

  end function GetQtyStuffForForwardModel
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ForwardModelVectorTools

! $Log$
! Revision 2.29  2018/02/27 01:33:46  livesey
! Fixed klutzy typo
!
! Revision 2.28  2018/02/27 00:51:01  livesey
! Added the supportedInstrumentModule functionality to the various search routines to support A-SMLS
!
! Revision 2.27  2014/07/18 23:14:29  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.26  2013/08/12 23:48:09  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.25  2013/06/12 02:20:02  vsnyder
! Cruft removal
!
! Revision 2.24  2013/04/03 23:23:33  vsnyder
! Don't look in config if it's not present
!
! Revision 2.23  2013/03/30 00:12:01  vsnyder
! Add GetQtyStuffForForwardModel
!
! Revision 2.22  2011/11/11 00:42:06  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.21  2011/06/16 20:19:47  vsnyder
! Add Frq argument to GetQuantityForForwardModel
!
! Revision 2.20  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.19  2009/04/20 16:36:26  pwagner
! Needed changes when identical types with different names allowed in L2PC files
!
! Revision 2.18  2008/10/03 16:26:47  livesey
! Added EXTINCTIONV2
!
! Revision 2.17  2008/09/29 22:56:54  vsnyder
! Add PRINT statement in Not_Used_Here to reduce compilation cascades
!
! Revision 2.16  2008/09/29 22:54:13  vsnyder
! Print radiometer name correctly in error message
!
! Revision 2.15  2006/07/17 20:06:37  livesey
! Added initialization of foundInFirst
!
! Revision 2.14  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.13  2004/11/01 20:19:35  vsnyder
! Moved QtyStuff_t and associated dump routine to ForwardModelConfig
!
! Revision 2.12  2004/10/16 17:28:28  livesey
! Added wasSpecific argument
!
! Revision 2.11  2004/10/07 23:25:50  vsnyder
! Move Dump_Qty_Stuff here from Get_Species_Data
!
! Revision 2.10  2004/06/10 00:59:56  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.9  2003/05/29 16:37:21  livesey
! New reflector argument to GetQuantityForForwardModel
!
! Revision 2.8  2003/05/05 23:00:24  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.7.2.1  2003/03/21 02:47:38  vsnyder
! Add a type with a pointer to a quantity, to make arrays of pointers
!
! Revision 2.7  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.6  2002/10/02 22:52:19  vsnyder
! Remove declaration for unused variable MATCHSC
!
! Revision 2.5  2002/10/02 22:51:07  vsnyder
! OOPS, mispelled 'parameter' as 'private'
!
! Revision 2.4  2002/10/02 22:49:13  vsnyder
! Cosmetic change to RCS stuff
!
! Revision 2.3  2002/09/26 21:41:23  livesey
! Moved error checking
!
! Revision 2.2  2002/09/26 18:01:34  livesey
! Bug fixes.
!
! Revision 2.1  2002/09/25 22:36:49  livesey
! First version
!
