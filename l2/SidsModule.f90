! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module SidsModule
  !=============================================================================

  ! This module evaluates the radiative transfer equation, and maybe
  ! its derivatives.  It is used for SIDS and L2PC runs.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: dump
  use ForwardModelConfig, only: ForwardModelConfig_T
  use ForwardModelWrappers, only: ForwardModel
  use ForwardModelIntermediate, only: ForwardModelIntermediate_T,&
    & ForwardModelStatus_T, DestroyForwardModelIntermediate
  use Init_Tables_Module, only: f_forwardModel, f_fwdModelIn, f_fwdModelExtra, &
    f_fwdModelOut, f_jacobian, f_perturbation
  use Lexer_Core, only: Print_Source
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MatrixModule_0, only: M_Absent, M_Full, M_Banded, M_Column_Sparse, &
    & MatrixElement_T
  use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
    GetFromMatrixDatabase, Matrix_Database_T, Matrix_T, DestroyBlock, CreateBlock, &
    & FindBlock
  use MoreTree, only: Get_Field_Id
  use Output_M, only: Output
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use VectorsModule, only: CopyVector, DestroyVectorInfo, Vector_T, operator(-)

  implicit none

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! subroutine SIDS ( Root, VectorDatabase, MatrixDatabase, FwdModelInfo )
  subroutine SIDS ( Root, VectorDatabase, MatrixDatabase, configDatabase)

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
    ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: configDatabase

    integer :: Error                    ! >= indicates an error occurred
    integer :: Field                    ! Of the "sids" specification
    integer :: config                   ! Index for config loop
    integer, dimension(:), pointer :: Configs ! Forward model configs
    type(vector_T), pointer :: FwdModelIn
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    type(vector_T), pointer :: Perturbation
    type(vector_T) :: SaveState         ! For numerical derivatives
    type(vector_T) :: Deviation         ! Result of perturbation
    type(vector_T) :: SaveResult        ! Result of forward model out
    integer :: I                        ! Subscript, loop inductor
    integer :: IxJacobian               ! Index of Jacobian in matrix database
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Son                      ! Of ROOT
    integer :: QUANTITY                 ! Index
    integer :: INSTANCE                 ! Index
    integer :: ELEMENT                  ! Index
    integer :: STATUS                   ! Flag
    logical :: DOTHISONE                ! Flag
    integer :: LOOPEND                  ! Loop limit
    integer :: COL                      ! Column in jacobian
    integer :: ROW                      ! Row in jacobian
    integer :: ROWINSTANCE              ! From jacobian
    integer :: ROWQUANTITY              ! From jacobian
    type (MatrixElement_T), pointer :: M0 ! A block from the jacobian

    type (ForwardModelIntermediate_T) :: ifm ! Work space for forward model
    type (ForwardModelStatus_T) :: fmStat ! Status for forward model

    real (r8) :: THISPTB                ! A scalar perturbation

    ! Error message codes
    integer, parameter :: NeedJacobian = 1   ! Needed if derivatives requested
    integer, parameter :: NotPlain = needJacobian + 1  ! Not a "plain" matrix
    integer, parameter :: PerturbationNotState = NotPlain + 1 ! Ptb. not same as state

    if ( toggle(gen) ) call trace_begin ( "SIDS", root )

    nullify ( configs, perturbation )

    ! Process the fields of the "sids" specification
    error = 0
    ixJacobian = 0

    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      select case ( field )
      case ( f_fwdModelIn )
        fwdModelIn => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_fwdModelExtra )
        fwdModelExtra => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_fwdModelOut )
        fwdModelOut => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_perturbation )
        perturbation => vectorDatabase(decoration(decoration(subtree(2,son))))
      case ( f_jacobian )
        ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
      case ( f_forwardModel )
        call Allocate_Test ( configs, nsons(son)-1, 'configs', ModuleName )
        do config = 2, nsons(son)
          configs(config-1) = decoration(decoration(subtree(config,son)))
        end do
      end select
    end do ! i = 2, nsons(root)

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
      endif
    end if

    ! Setup some stuff for the case where we're doing numerical derivatives.
    if ( associated ( perturbation ) ) then
      if ( perturbation%template%id /= fwdModelIn%template%id ) then
        call AnnounceError ( perturbationNotState )
        stop
      endif
      quantity = 1
      instance = 1
      element = 0
      loopEnd = perturbation%template%totalElements
    else
      ! Alternatively do simple one shot run.
      loopEnd = 1
      doThisOne = .true.
    endif

    ! Now have a possible loop over state vector elements, applying corrections.
    do i = 1, loopEnd

      ! Deal with the case where this is a perturbation.
      if ( associated ( perturbation ) ) then
        if ( element == 0 ) then
          ! Do a first unperturbed case
          call CopyVector ( saveState, fwdModelIn, clone=.true. )
          doThisOne=.true.
        else
          thisPtb =perturbation%quantities(quantity)%values(element,instance)
          if ( thisPtb /= 0.0 ) then
            call CopyVector ( fwdModelIn, saveState )
            fwdModelIn%quantities(quantity)%values(element,instance) = &
              & fwdModelIn%quantities(quantity)%values(element,instance) + &
              & thisPtb
            doThisOne=.true.
          else
            doThisOne=.false.
          endif
        end if
      end if

      ! Now loop over the MAFs / forward model configs and run the models
      if ( doThisOne ) then
        fmStat%newHydros = .true.
        fmStat%newScanHydros = .true.
        fmStat%maf = 0
        fmStat%finished = .false.
        if ( ixJacobian > 0 ) then
          call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', &
            & ModuleName )
        else ! because it's not optional in many places in the forward model
          call allocate_test ( fmStat%rows, 0, 'fmStat%rows', ModuleName )
        end if
        
        ! Loop over mafs
        do while ( .not. fmStat%finished )
          ! What if one config set finished but others still had more to do?
          ! Ermmm, think of this next time.
          fmStat%maf = fmStat%maf + 1
          do config = 1, size(configs)
            if ( ixJacobian > 0 .and. .not. associated(perturbation)) then
              call forwardModel ( configDatabase(configs(config)), &
                & FwdModelIn, FwdModelExtra, &
                & FwdModelOut, ifm, fmStat, Jacobian )
            else
              call forwardModel ( configDatabase(configs(config)), &
                & FwdModelIn, FwdModelExtra, &
                & FwdModelOut, ifm, fmStat )
            end if
          end do
        end do

        ! Tidy up the intermediate/status stuff from the model
        call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
        call DestroyForwardModelIntermediate ( ifm )
        
        ! Place the numerical derivative result into Jacobian if needed
        if ( associated ( perturbation ) ) then
          if ( element == 0 ) then
            call CopyVector (saveResult, fwdModelOut, clone=.true.)
          else
            deviation = fwdModelOut - saveResult
            ! Store the deviation in the jacobian
            col = FindBlock ( jacobian%col, quantity, instance )
            ! Loop over rows in the jacobian
            do row = 1, jacobian%row%nb
              rowQuantity = jacobian%row%quant(row)
              rowInstance = jacobian%row%inst(row)
              ! Is there anydeviaiton to store?
              if ( maxval ( abs ( &
                & deviation%quantities(rowQuantity)%values(:,rowInstance))) /= 0.0 ) then
                
                ! If so, this column of the block (creating if necessary)
                m0 => jacobian%block(row,col)
                if ( m0%kind == M_Absent ) &
                  & call CreateBlock ( jacobian, row, col, m_full )
                if ( any (m0%kind == (/ m_banded, m_column_sparse /) ) ) &
                  & call MLSMessage(MLSMSG_Error, ModuleName, &
                  & 'Unable to fill banded/column sparse blocks numerically')
                m0%values(:,element) = &
                  & deviation%quantities(rowQuantity)%values(:,rowInstance)/thisPtb

              end if                    ! Anything to place?
            end do                      ! Loop over rows
          endif                         ! Not very first run
        endif                           ! Doing perturbation.

      end if                            ! Do this element

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
    end if

    if ( toggle(gen) ) call trace_end ( "SIDS" )

    call deallocate_test ( configs, 'configs', ModuleName )

  contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code )
      integer, intent(in) :: Code       ! Index of error message

      error = max(error,1)
      call output ( '***** At ' )
      call print_source ( source_ref(root) )
      call output ( ' RetrievalModule complained: ' )
      select case ( code )
      case ( needJacobian )
        call output ( 'A Jacobian is required if derivatives are requested.', &
          & advance='yes' )
      case ( notPlain )
        call output ( 'The Jacobian matrix is not a "plain" matrix.', &
          & advance='yes' )
      case ( perturbationNotState )
        call output ( 'The perturbation vector is not the same type as fwdModelIn.', &
          & advance='yes' )
      end select
    end subroutine AnnounceError

  end subroutine SIDS

end module SidsModule

! $Log$
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
