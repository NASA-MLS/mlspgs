! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ForwardModelVectorTools          ! Tools for vectors in forward models

  ! This module contains routines needed to help a forward model get
  ! hold of the quantities it needs.

  implicit NONE

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  private

  public :: GetQuantityForForwardModel

contains

  ! ------------------------------ GetQuantityForForwardModel ----------------
  function GetQuantityForForwardModel ( vector, otherVector, quantityType, &
    & molecule, instrumentModule, radiometer, signal, sideband, &
    & molIndex, config, foundInFirst, noError )

    ! This function is in many senses like GetVectorQuantityByType, (to
    ! which it can revert), except that given a forwardModelConfig_T in
    ! config, and possibly an index into the molecules array, it can
    ! use the specificQuantities stuff in config to identify exactly
    ! the right quantity.
    use Intrinsic, only: L_VMR
    use Molecules, only: L_EXTINCTION
    use ForwardModelConfig, only: ForwardModelConfig_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: GetVectorQuantityByType, Vector_T, VectorValue_T
    use QuantityTemplates, only: QuantityTemplate_T
    use String_table, only: Get_String
    use MLSSignals_m, only: GetSignalName
    use MLSCommon, only: FindFirst
    use Intrinsic, only: Lit_Indices
    use Allocate_Deallocate, only: Allocate_test, Deallocate_test

    ! Dummy arguments
    type (Vector_T), target :: VECTOR ! First vector to look in
    type (Vector_T), target, optional :: OTHERVECTOR ! Second vector to look in
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in),  optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in),  optional :: INSTRUMENTMODULE ! Instrument module index
    integer, intent(in),  optional :: RADIOMETER   ! Radiometer index
    integer, intent(in),  optional :: SIGNAL       ! Signal index
    integer, intent(in),  optional :: SIDEBAND ! -1, 0, +1
    type (ForwardModelConfig_T), intent(in), optional :: CONFIG ! fwmConfig
    integer, intent(in),  optional :: MOLINDEX     ! Index into the molecules array
    logical, intent(out), optional :: FOUNDINFIRST ! Set if found in first vector
    logical, intent(in),  optional :: NOERROR ! Don't give error if not found
    ! Result
    type (VectorValue_T), pointer :: GetQuantityForForwardModel

    ! Local variables
    logical :: UseGetQuantityByType
    logical :: MyNoError                ! Copy of no error
    logical, dimension(:), pointer :: MATCHV ! Flags for vector 
    logical, dimension(:), pointer :: MATCHOV ! Flags for otherVector
    logical, dimension(:), pointer :: MATCHSC ! Flags for specificQuantities
    logical, dimension(:), pointer :: MATCH ! One of the match...s
    integer :: MolEntry                 ! How many similar molecules in list?
    integer :: NoFound                  ! Number of matches found so far.
    integer :: VectorIndex              ! Loop counter
    integer :: NoVectors                ! Number of vectors we've been given
    integer :: ThisMolecule             ! A molecule
    integer :: Quantity                 ! Loop counter
    type (Vector_T), pointer :: V       ! A vector
    type (QuantityTemplate_T), pointer :: QT ! A quantity template
    character(len=127) :: MSG

    ! Executable code

    ! Do a sanity check, get default values
    if ( present ( molIndex ) .and. .not. present ( config ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot have molIndex in GetQuantityForForwardModel without config" )
    myNoError = .false.
    if ( present ( noError ) ) myNoError = noError
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
        & quantityType, molecule, instrumentModule, radiometer, signal, sideband, &
        & foundInFirst, noError )
      return
    end if

    ! OK, looks like we have to do it the more complicated way.

    ! Setup a set of logicals to identify the matching quantities
    nullify ( matchV, matchOV )
    call Allocate_test ( matchV, size ( vector%quantities ), &
      & 'matchV', ModuleName )
    if ( present ( otherVector ) ) &
      & call Allocate_test ( matchOV, size ( otherVector%quantities ), &
      &   'matchOV', ModuleName )

    ! Now loop over the one or two vectors and work out which quantities might match.
    if ( present ( otherVector ) ) then
      noVectors = 2
    else
      noVectors = 1
    end if
    do vectorIndex = 1, noVectors
      if ( vectorIndex == 1 ) then
        v => vector
        match => matchV
      else
        v => otherVector
        match => matchOV
      end if
      do quantity = 1, size ( v%quantities )
        thisMolecule = 0
        qt => v%quantities(quantity)%template
        ! Now go through the quantities and see if they match
        match ( quantity ) = .false.
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
            if ( radiometer /= qt%radiometer .and. thisMolecule == l_extinction ) cycle
          else
            if ( radiometer /= qt%radiometer ) cycle
          end if
        end if
        if ( present(sideband) ) then
          if ( sideband /= qt%sideband ) cycle
        end if
        if ( present(signal) ) then
          if ( signal /= qt%signal ) cycle
        end if
        match ( quantity ) = .true.
      end do                            ! End loop over the quantities
    end do                              ! End loop over the one or two vectors

    ! Now decide exactly which ones we want
    if ( present ( molIndex ) ) then
      ! Some special thought for this case.
      ! Work out how many of the same molecules preceed or include this one.
      molEntry = count ( config%molecules ( 1:molIndex ) == &
        & config%molecules ( molIndex ) )
    else
      molEntry = 0
    end if

    ! First to see if any of the matches are on our 'specificQuantity'
    ! list, if so match them.  Though pay special attention in the molIndex case.
    noFound = 0
    if ( associated ( config%specificQuantities ) ) then
      specificVectorLoop: do vectorIndex = 1, noVectors
        if ( vectorIndex == 1 ) then
          v => vector
          match => matchV
        else
          v => otherVector
          match => matchOV
        end if
        do quantity = 1, size ( v%quantities )
          if ( match ( quantity ) .and. &
            & any ( v%template%quantities(quantity) == config%specificQuantities ) ) then
            noFound = noFound + 1
            if ( molEntry == 0 .or. molEntry == noFound ) then
              GetQuantityForForwardModel => v%quantities(quantity)
              if ( present ( foundInFirst ) ) foundInFirst = ( vectorIndex == 1 )
              exit specificVectorLoop
            end if
          end if
        end do
      end do specificVectorLoop
    end if

    if ( molEntry > 1 .and. noFound < molEntry ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unresolved ambiguity in molecule list' )

    ! Otherwise just get it from one or other vector
    if ( noFound == 0 ) then
      nonSpecificVectorLoop: do vectorIndex = 1, noVectors
        if ( vectorIndex == 1 ) then
          v => vector
          match => matchV
        else
          v => otherVector
          match => matchOV
        end if
        quantity = FindFirst ( match )
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
        call get_string ( lit_indices(radiometer), msg(len_trim(msg)+2:))
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

    call Deallocate_test ( matchV, 'matchV', ModuleName )
    call Deallocate_test ( matchOV, 'matchOV', ModuleName )

  end function GetQuantityForForwardModel

end module ForwardModelVectorTools

! $Log$
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
