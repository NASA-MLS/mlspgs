! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module SidsModule
  !=============================================================================

  ! This module evaluates the radiative transfer equation, and maybe
  ! its derivatives.  It is used for SIDS and L2PC runs.

  implicit none

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine SIDS ( Root, VectorDatabase, MatrixDatabase, HessianDatabase, &
    & configDatabase, chunk)

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Allocate
    use Chunks_M, only: MLSChunk_T
    use Expr_M, only: Expr
    use Forwardmodelconfig, only: Forwardmodelconfig_T
    use Forwardmodelwrappers, only: Forwardmodel
    use Forwardmodelintermediate, only: Forwardmodelstatus_T
    use Hessianmodule_1, only: Hessian_T, Inserthessianplane
    use HighOutput, only: OutputnamedValue
    use Init_Tables_Module, only: F_Destroyjacobian, F_Forwardmodel, &
      & F_Fwdmodelextra, F_Fwdmodelin, F_Fwdmodelout, &
      & F_Hessian, F_Jacobian, F_Mirrorhessian, &
      & F_Perturbation, F_SingleMAF, F_Switches, F_Tscat
    use, Intrinsic :: Iso_C_Binding, only: C_Intptr_T, C_Loc
    use Lexer_Core, only: Print_Source
    use MLSKinds, only: R8
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use Matrixmodule_0, only: M_Absent, M_Full, M_Banded, M_Column_Sparse, &
      & MatrixElement_T
    use Matrixmodule_1, only: Addtomatrix, Copymatrix, Createblock, &
      & Destroyblock, Destroymatrix, Dump, Findblock, &
      & GetfrommatrixDatabase, Matrix_Database_T, Matrix_T, Scalematrix
    use MLSL2options, only: L2Options
    use MLSL2timings, only: Add_To_Retrieval_Timing
    use MLSStringlists, only: Switchdetail
    use Moretree, only: Get_Field_Id, Get_Boolean
    use Output_M, only: Output
    use Scanmodelmodule, only: Destroyforwardmodelintermediate, &
      & Dumpinstancewindows
    use String_Table, only: Get_String
    use Time_M, only: Time_Now
    use Toggles, only: Gen, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Nsons, Subtree, Sub_Rosa, Where
    use Vectorsmodule, only: Vector_T, &
      & CopyVector, DestroyVectorinfo, Dump, Operator(-)

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
    ! Indexes an n_cf vertex
    type(vector_T), dimension(:), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(hessian_t), dimension(:), pointer :: HessianDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: configDatabase
    type(MLSChunk_T), intent(in) :: chunk

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Error                    ! >= indicates an error occurred
    integer :: ExprUnits(2)             ! From expr
    real(r8) :: ExprValue(2)            ! From expr
    integer :: Field                    ! Of the "sids" specification
    integer :: config                   ! Index for config loop
    integer, dimension(:), pointer :: Configs ! Forward model configs
    type(vector_T), pointer :: FwdModelIn
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    type(vector_T), pointer :: Perturbation
    type(vector_T) :: SaveState         ! For numerical derivatives
    type(matrix_T) :: SaveJacobian      ! The Jacobian matrix
    type(vector_T) :: Deviation         ! Result of perturbation
    type(vector_T) :: SaveResult        ! Result of forward model out
    integer :: I                        ! Subscript, loop inductor
    integer :: IxHessian                ! Index of Hessian in hessian database
    integer :: IxJacobian               ! Index of Jacobian in matrix database
    type(hessian_T), pointer :: Hessian ! The Hessian matrix
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Col                      ! Column in jacobian
    integer :: Element                  ! Index
    integer :: Instance                 ! Index
    integer :: LoopEnd                  ! Loop limit
    integer :: MAF1                     ! Loop limit
    integer :: MAF2                     ! Loop limit
    integer :: MAF                      ! Index
    integer :: Me = -1                  ! String index for trace
    integer :: Quantity                 ! Index
    logical :: GetAnaJac                ! Forward model needs to compute analytical Jacobians
    integer :: RowInstance              ! From jacobian
    integer :: RowQuantity              ! From jacobian
    integer :: Row                      ! Row in jacobian
    integer :: SingleMAF                ! From l2cf
    integer :: Son                      ! Of ROOT
    integer :: Status                   ! Flag
    integer :: SwitchLen                ! LEN_TRIM(Switches) on entry
    integer :: SwitchLenCur             ! LEN_TRIM(Switches) after command processing
    integer, allocatable :: ToKeep(:)   ! Molecule cross derivatives to keep
                                        ! after computing Hessian by
                                        ! perturbation. Union of molecule second
                                        ! derivatives from all configs.  Keep
                                        ! them all if this ends up empty as a
                                        ! result of no config specifying
                                        ! moleculeSecondDerivatives.
    logical :: DestroyJacobian          ! Flag
    logical :: DoThisOne                ! Flag
    logical :: DoTScat                  ! Flag
    logical :: MirrorHessian            ! Flag
    logical :: ShowPtb                  ! From -Sptb command-line option
    integer :: ShowPtbLevel             ! Number from -Sptb#
    real ::    T1
    type (MatrixElement_T), pointer :: M0 ! A block from the jacobian

    type (ForwardModelStatus_T) :: fmStat ! Status for forward model

    real (r8) :: THISPTB                ! A scalar perturbation

    ! Error message codes
    integer, parameter :: HessianNotJacobian = 1 ! Column vectors different
    integer, parameter :: NeedJacobian = HessianNotJacobian + 1 ! Needed if derivatives requested
    integer, parameter :: NotPlain = needJacobian + 1  ! Not a "plain" matrix
    integer, parameter :: PerturbationNotState = NotPlain + 1 ! Ptb. not same as state
    integer, parameter :: WrongDestroyJacobian = PerturbationNotState + 1 ! destroyJacobian with Hessian

    call trace_begin ( Me, "SIDS", root, cond=toggle(gen) )
    call time_now ( t1 )

    nullify ( configs, perturbation, Hessian )
    ! Process the fields of the "sids" specification
    doTScat = .false.
    error = 0
    ixJacobian = 0
    IxHessian = 0
    destroyJacobian = .false.
    mirrorHessian = .false.
    singleMAF = -1
    fwdModelExtra => NULL()             ! Can be omitted
    switchLen = len_trim(switches)
    switchLenCur = switchLen + 1
    switches(switchLenCur:switchLenCur) = ','

    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      select case ( field )
      case ( f_destroyJacobian )
        destroyJacobian = Get_boolean(son)
      case ( f_mirrorHessian )
        mirrorHessian = Get_boolean(son)
      case ( f_forwardModel )
        call Allocate_Test ( configs, nsons(son)-1, 'configs', ModuleName )
        do config = 2, nsons(son)
          configs(config-1) = decoration(decoration(subtree(config,son)))
        end do
      case ( f_fwdModelExtra )
        fwdModelExtra => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_fwdModelIn )
        fwdModelIn => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_fwdModelOut )
        fwdModelOut => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_jacobian )
        ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
      case ( f_hessian )
        ixHessian =  decoration(subtree(2,son)) ! hessian: hessian vertex
      case ( f_perturbation )
        perturbation => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_singleMAF )
        call expr ( subtree(2,son), exprUnits, exprValue )
        singleMAF = exprValue(1)
      case ( f_switches )
        call get_string ( sub_rosa(subtree(2,son)), switches(switchLenCur+1:), strip=.true. )
        switchLenCur = len_trim(switches) + 1
        switches(switchLenCur:switchLenCur) = ','
      case ( f_TScat )
        doTScat = Get_boolean(son)
      end select
    end do ! i = 2, nsons(root)

    showPtbLevel = switchDetail( switches, 'ptb' )
    showptb = showPtbLevel > -1 .and. associated( perturbation )

    ! Now if we weren't given a forward model extra, point it to the forwardModelIn
    ! This seems rather cheesy.  However, the alternative is changing a lot of
    ! code from having FwdModelExtra as optional to pointer.
    if ( .not. associated ( fwdModelExtra ) ) &
      & fwdModelExtra => fwdModelIn

    ! Check that we have a Jacobian if we need one
    if ( ixJacobian > 0 ) then
      i = decoration(ixJacobian)
      call getFromMatrixDatabase ( matrixDatabase(i), jacobian )
      if ( .not. associated(jacobian) ) call announceError ( notPlain )
    else
      if ( any( (/configDatabase(configs)%atmos_Der, &
        &         configDatabase(configs)%spect_Der, &
        &         configDatabase(configs)%temp_der/) ) ) then
        call announceError ( needJacobian )
      end if
    end if

    ! Check we have a Jacobian if we have a Hessian
    if ( ixHessian > 0 .and. ixJacobian <= 0 ) call announceError( needJacobian )
    if ( ixHessian > 0 ) then
      hessian => HessianDatabase ( -decoration ( ixHessian ) )
      if ( hessian%col%vec%template%name /= jacobian%col%vec%template%name &
        & .or. hessian%row%vec%template%name /= jacobian%row%vec%template%name &
        & .or. (hessian%col%instFirst .neqv. jacobian%col%instFirst) &
        & .or. (hessian%row%instFirst .neqv. jacobian%row%instFirst) )&
        & call announceError ( HessianNotJacobian )
    end if

    ! Check that if we set destroyJacobian we're not asking for Hessian information
    if ( destroyJacobian .and. ( ixHessian > 0 ) ) call announceError ( wrongDestroyJacobian )

    ! Now work out if we're going to be asking the forward model to compute analytical Jacobians
    ! or if it's up to us to do so
    getAnaJac = ( ixJacobian > 0 ) .and. ( .not. associated ( perturbation ) .or. &
      & ( ixHessian > 0 ) )

    ! Setup some stuff for the case where we're doing numerical derivatives.
    if ( associated ( perturbation ) ) then
      if ( perturbation%template%name /= fwdModelIn%template%name ) &
        & call AnnounceError ( perturbationNotState )
      quantity = 1
      instance = 1
      element = 0
      loopEnd = perturbation%template%totalElements
      if ( showptb ) then
        call outputnamedValue( 'loopEnd', loopEnd )
        call outputnamedValue( 'size(configs)', size(configs) )
      end if
    else
      ! Alternatively do simple one shot run.
      loopEnd = 1
      doThisOne = .true.
    end if

    if ( error > 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Errors prevent SIDS run.' )

    configDatabase(configs)%generateTScat = doTScat

    if ( ixHessian > 0 ) call get_toKeep ! Cross derivatives to keep in Hessian

    ! Now have a possible loop over state vector elements, applying corrections.
    do i = 1, loopEnd

      ! Deal with the case where this is a perturbation.
      if ( associated ( perturbation ) ) then
        if ( element == 0 ) then
          ! Do a first unperturbed case.  Save the initial state first.
          call CopyVector ( saveState, fwdModelIn, clone=.true. )
          doThisOne=.true.
        else
          thisPtb =perturbation%quantities(quantity)%values(element,instance)
          doThisOne = thisPtb /= 0.0
          if ( doThisOne ) then
            call CopyVector ( fwdModelIn, saveState ) ! Get saved state
            ! Perturb it
            fwdModelIn%quantities(quantity)%values(element,instance) = &
              & fwdModelIn%quantities(quantity)%values(element,instance) + &
              & thisPtb
          end if
        end if
      end if

      ! Now loop over the MAFs / forward model configs and run the models
      if ( doThisOne ) then
        fmStat%newScanHydros = .true.

        if ( ixJacobian > 0 ) then
          call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
            & ModuleName )
        else ! because it's not optional in many places in the forward model
          call allocate_test ( fmStat%rows, 0, 'fmStat%rows', ModuleName )
        end if

        if ( showptb .and. instance == 1 ) then
          call dump( perturbation%quantities(quantity), details=showPtbLevel-1, &
            & name='Perturbation', vector=perturbation )
          call dump ( saveState%quantities(quantity), details=showPtbLevel-1, &
            & name='Save State', vector=saveState )
        end if

        ! Loop over forward model configs
        do config = 1, size(configs)
          ! Work out the loop limits
          if ( singleMAF == -1 ) then
            maf1 = 1
            maf2 = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
            if (configDatabase(configs(config))%skipOverlaps) then
              maf1 = maf1 + chunk%noMAFsLowerOverlap
              maf2 = maf2 - chunk%noMAFsUpperOverlap
            end if
          else
            if ( (singleMAF < 1) .or. &
              & (singleMAF > chunk%lastMAFIndex - chunk%firstMAFIndex + 1) ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Illegal value of singleMAF' )
            maf1 = singleMAF
            maf2 = singleMAF
          end if
          if ( switchDetail( switches, 'fiw' ) > -1 ) &
            & call DumpInstanceWindows ( &
            & fwdModelIn, fwdModelExtra, maf1, maf2, &
            & configDatabase(configs(config)), &
            & details=switchDetail( switches, 'fiw' ) &
            & )
          if ( L2Options%SkipRetrieval ) cycle
          ! Loop over mafs
          do maf = maf1, maf2
            fmStat%maf = maf
              call add_to_retrieval_timing( 'sids', t1 )
            if ( getAnaJac ) then
              if ( associated ( Hessian ) .and. .not. associated(perturbation) ) then
                call forwardModel ( configDatabase(configs(config)), &
                  & FwdModelIn, FwdModelExtra, &
                  & FwdModelOut, fmStat, Jacobian, Hessian, vectorDatabase )
              else
                call forwardModel ( configDatabase(configs(config)), &
                  & FwdModelIn, FwdModelExtra, &
                  & FwdModelOut, fmStat, Jacobian, vectors=vectorDatabase )
              end if
            else
              call forwardModel ( configDatabase(configs(config)), &
                & FwdModelIn, FwdModelExtra, &
                & FwdModelOut, fmStat, vectors=vectorDatabase )
            end if
            call time_now ( t1 )
          end do                        ! MAF loop

          ! Destroy jacobian if asked to
          if (destroyJacobian .and. ixJacobian > 0 ) then
            call DestroyBlock ( Jacobian )
            allocate ( Jacobian%block ( jacobian%row%nb, jacobian%col%nb ), &
              & STAT=status )
            addr = 0
            if ( status == 0 ) then
              if ( size(Jacobian%block) > 0 ) &
                & addr = transfer(c_loc(Jacobian%block(1,1)), addr)
            end if
            call test_allocate ( status, moduleName, 'jacobian%block', &
              & uBounds=[jacobian%row%nb,jacobian%col%nb], &
              & elementSize=storage_size(Jacobian%block)/8, address=addr )
          end if
        end do                          ! Forward model config loop

        ! Tidy up the intermediate/status stuff from the model
        call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
        call DestroyForwardModelIntermediate ! in case scan model got used
        fmStat%newScanHydros = .true.
        if ( L2Options%SkipRetrieval ) cycle

        ! Place the numerical derivative result into Jacobian if needed
        if ( associated ( perturbation ) ) then
          ! Where to store deviation in Jacobian or Hessian
          col = FindBlock ( jacobian%col, quantity, instance )
          if ( element == 0 ) then
            call CopyVector (saveResult, fwdModelOut, clone=.true.)
            if ( ixHessian > 0 ) call CopyMatrix ( saveJacobian, jacobian )
          else
            if ( ixHessian <= 0 ) then
              ! We're doing perturbations to get Jacobians
              deviation = fwdModelOut - saveResult
              ! Loop over rows in the jacobian
              do row = 1, jacobian%row%nb
                rowQuantity = jacobian%row%quant(row)
                rowInstance = jacobian%row%inst(row)
                ! Is there any deviation to store?
                if ( maxval ( abs ( &
                   & deviation%quantities(rowQuantity)%values(:,rowInstance))) /= 0.0 ) then

                  ! If so, this column of the block (creating if necessary)
                  m0 => jacobian%block(row,col)
                  if ( m0%kind == M_Absent ) then
                    call CreateBlock ( jacobian, row, col, m_full )
                    m0%values = 0.0
                  end if
                  if ( any (m0%kind == (/ m_banded, m_column_sparse /) ) ) &
                    & call MLSMessage(MLSMSG_Error, ModuleName, &
                    & 'Unable to fill banded/column sparse blocks numerically')
                  m0%values(:,element) = &
                    & deviation%quantities(rowQuantity)%values(:,rowInstance)/thisPtb

                end if                    ! Anything to place?
              end do                      ! Loop over rows
            else
              ! We're doing perturbations to get Hessians
              call AddToMatrix ( jacobian, saveJacobian, -1.0_r8 )
              call ScaleMatrix ( jacobian, 1.0/thisPtb )
              ! Insert this difference
              if ( switchDetail( switches, 'hess' ) > 0 ) then
                call output ( 'Inserting Jacobian numerical deriv. into Hessian', advance='yes' )
                call dump( jacobian, details=1 )
              end if
              call InsertHessianPlane ( hessian, jacobian, col, element, &
                & toKeep, mirror=mirrorHessian )
            end if
          end if                         ! Not very first run
        end if                           ! Doing perturbation.

      end if                             ! Do this element

      ! Work out which element is next if we're perturbing
      if ( associated ( perturbation ) ) then
        element = element + 1
        if ( element == &
          & perturbation%quantities(quantity)%template%instanceLen + 1) then
          element = 1
          instance = instance + 1
          if ( instance == &
            & perturbation%quantities(quantity)%template%noInstances + 1 ) then
            instance = 1
            quantity = quantity + 1
          end if
        end if
      end if

    end do                              ! Possible loop over perturbation elements

    if ( associated ( perturbation ) ) then
      call CopyVector ( fwdModelIn, saveState )
      call CopyVector ( fwdModelOut, saveResult )
      call DestroyVectorInfo ( saveState )
      call DestroyVectorInfo ( saveResult )
      call DestroyVectorInfo ( deviation )
      if ( ixHessian > 0 ) then
        call CopyMatrix ( jacobian, saveJacobian )
        call DestroyMatrix ( saveJacobian )
      end if
    end if

    configDatabase(configs)%generateTScat = .false.

    call deallocate_test ( configs, 'configs', ModuleName )
    call add_to_retrieval_timing( 'sids', t1 )

    switches(switchLen+1:) = '' ! Clobber switches from SIDS command

    call trace_end ( "SIDS", cond=toggle(gen) )

  contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code )
      integer, intent(in) :: Code       ! Index of error message

      error = max(error,1)
      call output ( '***** At ' )
      call print_source ( where(root) )
      call output ( ' SidsModule complained: ' )
      select case ( code )
      case ( HessianNotJacobian )
        call output ( 'Hessian and Jacobian are not compatible.', advance='yes' )
      case ( needJacobian )
        call output ( 'A Jacobian is required if derivatives are requested.', &
          & advance='yes' )
      case ( notPlain )
        call output ( 'The Jacobian matrix is not a "plain" matrix.', &
          & advance='yes' )
      case ( perturbationNotState )
        call output ( 'The perturbation vector is not the same type as fwdModelIn.', &
          & advance='yes' )
      case ( wrongDestroyJacobian )
        call output ( 'Cannot set destroyJacobian and ask for a Hessian', advance='yes' )
      end select
    end subroutine AnnounceError

    ! -----------------------------------------------  Get_ToKeep  -----
    subroutine Get_ToKeep
      integer :: Config, I ! Loop index
      integer :: ToKeepGuess(sum( (/ ( count(configDatabase(configs(config))%moleculeSecondDerivatives), &
                                     & config = 1, size(configs) ) /) ) )
      i = 0
      do config = 1, size(configs)
        toKeepGuess(i+1:i+count(configDatabase(configs(config))%moleculeSecondDerivatives)) = &
          & pack(configDatabase(configs(config))%molecules, &
          &      configDatabase(configs(config))%moleculeSecondDerivatives)
        i = i + count(configDatabase(configs(config))%moleculeSecondDerivatives)
      end do
      ! Delete duplicates
      do i = 2, size(toKeepGuess)
        if ( any( toKeepGuess(1:i-1) == toKeepGuess(i) ) ) toKeepGuess(i) = 0
      end do
      allocate ( toKeep(count(toKeepGuess>0)), stat=status )
      call test_allocate ( status, moduleName, 'ToKeep' )
      toKeep = pack(toKeepGuess,toKeepGuess>0)
    end subroutine Get_ToKeep

  end subroutine SIDS

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SidsModule

! $Log$
! Revision 2.78  2018/09/13 20:24:50  pwagner
! Moved changeable options to new L2Options; added DumpOptions
!
! Revision 2.77  2017/12/07 01:01:24  vsnyder
! Don't use host-associated variable as a DO index
!
! Revision 2.76  2017/07/27 01:41:50  vsnyder
! Better debug printing
!
! Revision 2.75  2015/03/28 02:51:19  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.74  2014/09/05 01:20:54  vsnyder
! Remove USE for unreferenced name
!
! Revision 2.73  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.72  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.71  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.70  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.69  2013/08/30 02:45:48  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.68  2013/03/15 20:36:03  vsnyder
! Add 'switches' field to SIDS
!
! Revision 2.67  2013/03/01 01:09:50  pwagner
! 'fiw' switch dumps instance window
!
! Revision 2.66  2012/10/11 21:02:16  pwagner
! Fix timing error--dont include ForwardModel times
!
! Revision 2.65  2011/05/17 00:44:14  vsnyder
! Remove ill-advised handling of perturbation for logBasis molecules.  The
! logBasis flag only affects copying to the forward model's grid structures,
! not how mixing ratios are stored in vector quantities.
!
! Revision 2.64  2011/05/17 00:26:25  vsnyder
! Handle log-basis perturbation correctly
!
! Revision 2.63  2011/04/02 01:23:22  vsnyder
! Make a list of the molecule cross derivatives to keep from the union of
! all molecules specified in moleculeSecondDerivatives fields of all configs.
!
! Revision 2.62  2011/03/30 00:41:07  vsnyder
! Don't ask for analytic Hessian if using perturbation
!
! Revision 2.61  2011/03/25 20:44:45  vsnyder
! Initially nullify Hessian.  Use Test_Allocate.  Delete cross derivatives
! that are not requested, when storing a Hessian plane.
!
! Revision 2.60  2010/09/17 00:08:18  pwagner
! Should not crash if 'ptb' switch set but no perturbation asoociated
!
! Revision 2.59  2010/08/27 06:36:41  yanovsky
! Can now call forwardModel with Hessian as argument
!
! Revision 2.58  2010/08/06 23:04:24  pwagner
! Added 'hess' switch
!
! Revision 2.57  2010/05/13 23:48:19  pwagner
! Added -Sptb to show perturbation
!
! Revision 2.56  2010/03/29 18:40:09  pwagner
! Repaired error due to undefined ixHessian
!
! Revision 2.55  2010/03/24 20:56:11  vsnyder
! Minor tweaks in Hessian generation
!
! Revision 2.54  2010/02/25 18:20:12  pwagner
! Adds support for new Hessian data type
!
! Revision 2.53  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.52  2008/06/06 01:59:28  vsnyder
! Set DoTScat flag in configs, some cannonball polishing
!
! Revision 2.51  2007/06/29 19:32:07  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.50  2007/04/03 17:39:36  vsnyder
! Replace pointer attribute on VectorDatabase with target attribute
!
! Revision 2.49  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.48  2004/12/27 23:06:15  vsnyder
! Remove unreferenced use names
!
! Revision 2.47  2004/06/21 22:50:55  livesey
! Some trapping of silly values of singleMAF in sids command
!
! Revision 2.46  2004/05/19 19:16:12  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.45  2004/03/31 03:59:43  livesey
! Added singleMAF option
!
! Revision 2.44  2003/09/11 23:16:36  livesey
! Now hands the vectors database to the forward model to support the
! linearized forward model's xStar / yStar capabilities.
!
! Revision 2.43  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.42  2002/10/08 17:36:23  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.41  2002/09/23 18:04:47  livesey
! Got rid of unnecessary timing
!
! Revision 2.40  2002/09/18 23:56:01  vsnyder
! Call time_now at end of add_to_retrieval_timing
!
! Revision 2.39  2002/01/09 00:51:05  livesey
! Allow fwdModelExtra to be omitted, a different way
!
! Revision 2.38  2002/01/08 18:15:57  livesey
! Made fwdModelExtra optional
!
! Revision 2.37  2001/11/27 23:34:49  pwagner
! Split forward model timings into four types
!
! Revision 2.36  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.35  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models
!
! Revision 2.34  2001/10/01 22:54:22  pwagner
! Added subsection timings for Retrieval section
!
! Revision 2.33  2001/05/26 00:21:38  livesey
! Added zeroing of jacobian blocks in numerical derivative case.
!
! Revision 2.32  2001/05/10 01:08:11  livesey
! Added destroyJacobian option
!
! Revision 2.31  2001/05/09 01:50:31  vsnyder
! Make sure fmstat%rows is defined (even though it's not used)
!
! Revision 2.30  2001/05/08 21:33:56  livesey
! Removed some print statements
!
! Revision 2.29  2001/05/05 00:02:27  livesey
! Got a numerical derivative code working.
!
! Revision 2.28  2001/05/03 23:09:09  livesey
! Removed temporary code to delete jacobian
!
! Revision 2.27  2001/05/03 20:33:51  vsnyder
! Added a nullify and did some cosmetic changes
!
! Revision 2.26  2001/05/01 23:52:54  vsnyder
! Allocate and deallocate fmStat%rows here
!
! Revision 2.25  2001/05/01 06:57:30  livesey
! *** empty log message ***
!
! Revision 2.24  2001/05/01 00:20:34  livesey
! Sets up fmStat%rows correctly.
!
! Revision 2.23  2001/04/26 19:48:11  livesey
! Now uses ForwardModelWrappers
!
! Revision 2.22  2001/04/26 00:57:20  vsnyder
! Deallocate fmStat%rows, cosmetic changes
!
! Revision 2.21  2001/04/24 23:11:40  vsnyder
! Remove 'Done forward model!' print
!
! Revision 2.20  2001/04/19 23:56:01  livesey
! New fmStat
!
! Revision 2.19  2001/04/19 20:30:24  livesey
! Added call to DestroyForwardModelIntermediate
!
! Revision 2.18  2001/04/12 21:42:24  livesey
! Another interim version, forgot to nullify a pointer.
!
! Revision 2.17  2001/04/12 18:13:28  vsnyder
! OOPS! Hadn't saved it from the editor!
!
! Revision 2.16  2001/04/12 17:48:11  livesey
! Moved maf increment in from ForwardModel to here.
!
! Revision 2.15  2001/04/10 23:28:10  livesey
! Interim working version
!
! Revision 2.14  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.13  2001/04/10 00:24:30  vsnyder
! Add an error message if Jacobian isn't 'plain'
!
! Revision 2.12  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.11  2001/03/30 00:07:24  livesey
! Removed FMC in call to forwardModel
!
! Revision 2.10  2001/03/28 23:47:33  livesey
! Made it so it can run in a loop
!
! Revision 2.9  2001/03/25 00:50:31  livesey
! Interim version, bug with frequency averaging
!
! Revision 2.8  2001/03/20 02:30:15  livesey
! Interim version, gets same numbers as Zvi
!
! Revision 2.7  2001/03/17 00:45:28  livesey
! Moved to new ForwardModelConfig_T
!
! Revision 2.6  2001/03/09 01:35:18  vsnyder
! Don't run fwd model if derivatives requested but Jacobian not supplied
!
! Revision 2.5  2001/03/09 01:30:10  vsnyder
! Don't ask for Jacobian if convolution is requested
!
! Revision 2.4  2001/03/08 23:52:24  vsnyder
! Process sons correctly
!
! Revision 2.3  2001/03/08 20:11:19  zvi
! *** empty log message ***
!
! Revision 2.2  2001/03/08 03:23:09  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.1  2001/03/08 00:00:08  vsnyder
! Initial commit.
!
