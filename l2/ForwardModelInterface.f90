! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
!=============================================================================

! Set up the forward model.  Interface from the retrieve step to the
! forward model.

!??? Do we want a forward model database ???

  use Expr_M, only: EXPR
  use Init_Tables_Module, only: field_first, field_last
  use Init_Tables_Module, only: F_ATMOS_DER, F_DO_CONV, F_DO_FREQ_AVG, &
    & F_FREQUENCY, F_SPECT_DER, F_TEMP_DER
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Node_ID, Nsons, Source_Ref, Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: Vector_T

  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO

  implicit NONE
  private
  public :: ForwardModel, ForwardModelGlobalSetup, ForwardModelInfo_T, &
    & ForwardModelSetup

  type ForwardModelInfo_T
    logical :: Atmos_Der      ! Do atmospheric derivatives
    logical :: Do_Conv        ! Do convolution
    logical :: Do_Freq_Avg    ! Do Frequency averaging
    real(r8) :: The_Freq      ! Frequency to use if .not. do_freq_avg
    logical :: Spect_Der      ! Do spectroscopy derivatives
    logical :: Temp_Der       ! Do temperature derivatives
  end type ForwardModelInfo_T

  integer :: Error            ! Error level -- 0 = OK

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root, ForwardModelInfo )
  ! Process the forwardModel specification to produce ForwardModelInfo.

    integer :: Root                     ! of the forwardModel specification.
                                        ! Indexes either a "named" or
                                        ! "spec_args" vertex.
    type(forwardModelInfo_T), intent(inout) :: ForwardModelInfo

    integer :: Field                    ! Field index -- f_something
    integer :: I                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    integer :: Son                      ! Some subtree of root.
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ForwardModelGlobalSetup", root )
    if ( node_id(root) == n_named ) then
      name = subtree(1, root)
      key = subtree(2, root)
    else
      name = 0
      key = root
    end if

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      select case ( field )
      case ( f_atmos_der )
        forwardModelInfo%atmos_der = get_boolean(son)
      case ( f_do_conv )
        forwardModelInfo%do_conv = get_boolean(son)
      case ( f_do_freq_avg )
        forwardModelInfo%do_freq_avg = get_boolean(son)
      case ( f_frequency )
        call expr ( subtree(2,son), units, value, type )
        forwardModelInfo%the_freq = value(1)
      case ( f_spect_der )
        forwardModelInfo%spect_der = get_boolean(son)
      case ( f_temp_der )
        forwardModelInfo%temp_der = get_boolean(son)
      case default
        ! Shouldn't get here if the type checker worked
      end select
    end do ! i = 2, nsons(key)
    if ( toggle(gen) ) call trace_end ( "ForwardModelGlobalSetup" )
  end subroutine ForwardModelGlobalSetup

  ! ------------------------------------------  ForwardModelSetup  -----
  subroutine ForwardModelSetup ( Root, VectorDatabase, MatrixDatabase, &
    &                            ForwardModelInfo )
  ! Process the forwardModel specification to produce ForwardModelInfo.

    integer :: Root                     ! of the forwardModel specification.
                                        ! Indexes either a "named" or
                                        ! "spec_args" vertex.
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelInfo_T), intent(out) :: ForwardModelInfo

    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    integer :: Son                      ! Some subtree of root.

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ForwardModelSetup", root )
    if ( node_id(root) == n_named ) then
      name = subtree(1, root)
      key = subtree(2, root)
    else
      name = 0
      key = root
    end if

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    got = .false.
    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case default
        ! Shouldn't get here if the type checker worked
      end select
    end do ! i = 2, nsons(key)
    if ( toggle(gen) ) call trace_end ( "ForwardModelSetup" )

  end subroutine ForwardModelSetup

  ! -----------------------------------------------  ForwardModel  -----
! subroutine ForwardModel ( FwdModelInfo, FwdModelExtra, FwdModelIn, &
!   &                       Jacobian, RowBlock, FwdModelOut )
  subroutine ForwardModel ( FwdModelInfo, FwdModelExtra, FwdModelIn, &
    &                       Jacobian, RowBlock, FwdModelOut, FMC, FMI, TFMI )
    type(forwardModelInfo_T), intent(in) :: FwdModelInfo ! From ForwardModelSetup
    type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn ! ???
    type(matrix_T), intent(inout), optional :: Jacobian
    integer, intent(in), optional :: RowBlock          ! With which block of
    ! rows of F and Jacobian are we computing? All of them if absent.
    type(vector_T), intent(inout), optional :: FwdModelOut  ! Radiances, etc.

!??? Begin temporary stuff to start up the forward model
  type(fwd_mdl_config), optional :: FMC
  type(fwd_mdl_info), dimension(:), pointer, optional :: FMI
  type(temporary_fwd_mdl_info), dimension(:), pointer, optional :: TFMI
!??? End of temporary stuff to start up the forward model

  end subroutine ForwardModel

! =====     Private Procedures     =====================================
  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( Code, Where, FieldIndex )
    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: Where      ! Where in the tree did the error occur?
    integer, intent(in) :: FieldIndex ! f_...

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref ( where ) )
    call output ( ' ForwardModelSetup complained: ' )
    select case ( code )
    end select
  end subroutine AnnounceError
end module ForwardModelInterface

! $Log$
! Revision 2.5  2001/03/08 00:42:09  vsnyder
! Add temporary stuff to use with L2_Load
!
! Revision 2.4  2001/03/07 23:59:52  vsnyder
! Add stuff for SIDS.
!
! Revision 2.3  2001/02/21 00:07:57  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.2  2001/02/08 00:56:11  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.1  2001/02/07 00:52:27  vsnyder
! Initial commit
!
