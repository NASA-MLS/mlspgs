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
  use Dump_0, only: dump
  use ForwardModelIntermediate, only: ForwardModelIntermediate_T,&
    & ForwardModelStatus_T
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
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
    integer, pointer, dimension(:) :: configs ! Forward model configs
    type(vector_T), pointer :: FwdModelIn
    type(vector_T), pointer :: FwdModelExtra
    type(vector_T), pointer :: FwdModelOut
    integer :: I                        ! Subscript, loop inductor
    integer :: IxJacobian               ! Index of Jacobian in matrix database
    type(matrix_T), pointer :: Jacobian ! The Jacobian matrix
    integer :: Son                      ! Of ROOT

    type (ForwardModelIntermediate_T) :: ifm ! Work space for forward model
    type (ForwardModelStatus_T) :: fmStat ! Status for forward model

    ! Error message codes
    integer, parameter :: NeedJacobian = 1   ! Needed if derivatives requested
    integer, parameter :: NotPlain = needJacobian + 1  ! Not a "plain" matrix

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
        call Allocate_Test(configs, nsons(son)-1, 'configs', ModuleName)
        do config = 1, nsons(son)-1
          configs(config) = decoration(decoration(subtree(config+1,son)))
        end do
      end select
    end do ! i = 2, nsons(root)

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
    
    fmStat%newHydros = .true.
    fmStat%maf = 1

    ! Loop over mafs
    do
      do config = 1, size(configs)
        if ( ixJacobian > 0 ) then
          call forwardModel ( configDatabase(configs(config)), &
            & FwdModelIn, FwdModelExtra, &
            & FwdModelOut, ifm, fmStat, Jacobian)
        else
          call forwardModel ( configDatabase(configs(config)), &
            & FwdModelIn, FwdModelExtra, &
            & FwdModelOut, ifm, fmStat)
        end if
      end do
      ! What if one config set finished but others still had more to do?
      ! Ermmm, think of this next time.
      if (fmStat%finished) exit
      fmStat%maf = fmStat%maf + 1
    end do

    if ( toggle(gen) ) call trace_end ( "SIDS" )

    call deallocate_test( configs, 'configs', ModuleName)

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
      end select
    end subroutine AnnounceError

  end subroutine SIDS

end module SidsModule

! $Log$
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
