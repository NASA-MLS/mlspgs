! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module SidsModule
  !=============================================================================

  ! This module evaluates the radiative transfer equation, and maybe
  ! its derivatives.  It is used for SIDS and L2PC runs.

  use ForwardModelConfig, only: ForwardModelConfig_T
  use ForwardModelInterface, only: ForwardModel
  use Init_Tables_Module, only: f_forwardModel, f_fwdModelIn, f_fwdModelExtra, &
    f_fwdModelOut, f_jacobian
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, &
    GetFromMatrixDatabase, Matrix_Database_T, Matrix_T
  use MoreTree, only: Get_Field_Id
  use Output_M, only: Output
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use VectorsModule, only: Vector_T
  use l2_load_m, only: l2_load
  use Dump_0, only: dump

  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! subroutine SIDS ( Root, VectorDatabase, MatrixDatabase, FwdModelInfo )
  subroutine SIDS ( Root, VectorDatabase, MatrixDatabase, configDatabase, &
    & FMC, FMI, TFMI )

    ! Dummy arguments:
    integer, intent(in) :: Root         ! Of the relevant subtree of the AST
    ! Indexes an n_cf vertex
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelConfig_T), dimension(:), pointer :: configDatabase

    !??? Begin temporary stuff to start up the forward model
    type(fwd_mdl_config) :: FMC
    type(fwd_mdl_info), dimension(:), pointer :: FMI
    type(temporary_fwd_mdl_info), dimension(:), pointer :: TFMI
    !??? End of temporary stuff to start up the forward model

    type (ForwardModelConfig_T), pointer :: CONFIG ! Selected configuration
    integer :: Error                    ! >= indicates an error occurred
    integer :: Field                    ! Of the "sids" specification
    type(vector_T), pointer :: FwdModelIn
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    integer :: I                        ! Subscript, loop inductor
    integer :: IxJacobian               ! Index of Jacobian in matrix database
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Son                      ! Of ROOT

    ! Error message codes
    integer, parameter :: NeedJacobian = 1   ! Needed if derivatives requested

    if ( toggle(gen) ) call trace_begin ( "SIDS", root )

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
      case ( f_jacobian )
        ixJacobian = decoration(subtree(2,son)) ! jacobian: matrix vertex
      case ( f_forwardModel )
        config => configDatabase(decoration(decoration(subtree(2,son))))
      end select
    end do ! i = 2, nsons(root)

    if ( ixJacobian > 0 ) then
      i = decoration(ixJacobian)
      call getFromMatrixDatabase ( matrixDatabase(i), jacobian )
      call forwardModel ( config, FwdModelExtra, FwdModelIn, &
        &                   Jacobian, FwdModelOut=FwdModelOut, &
        &                   FMI=FMI(1), TFMI=TFMI(1)) !???  temporary
      !     &                  FMI=FMI,TFMI=TFMI) !??? Last line temporary
    else if ( config%atmos_Der .or. config%spect_Der .or. config%temp_der ) then
      call announceError ( needJacobian )
    else
      print*,'Calling forward model without derivatives'
      call forwardModel ( config, FwdModelExtra, FwdModelIn, &
        &                   FwdModelOut=FwdModelOut, &
        &                   FMI=FMI(1), TFMI=TFMI(1)) !??? temporary
      !     &               FMI=FMI,TFMI=TFMI) !??? Last line temporary
      print*,'Got back from forward model'
    end if
    print*,'Done the forward model!!!'
    if ( toggle(gen) ) call trace_end ( "SIDS" )

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
      end select
    end subroutine AnnounceError

  end subroutine SIDS

end module SidsModule

! $Log$
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
