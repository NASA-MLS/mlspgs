module GLOBAL_SETTINGS

  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT, &
    & S_FORWARDMODEL
  use MLSCommon, only: R8
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE

!??? Begin temporary stuff to start up the forward model
  use ForwardModelInterface, only: AddForwardModelConfigToDatabase, &
    ConstructForwardModelConfig, ForwardModelGlobalSetup, ForwardModelConfig_T
  use INIT_TABLES_MODULE, only: F_ZVI, S_ForwardModelGlobal, S_L2LOAD
  use L2_Load_M, only: L2_Load
  use L2_test_structures_m, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use String_Table, only: Get_String
  use TREE, only: NODE_ID, DECORATE
  use TREE_TYPES, only: N_EQUAL, N_NAMED
!??? End of temporary stuff to start up the forward model

  implicit NONE

  private

  public :: SET_GLOBAL_SETTINGS

  logical, public :: ALLOW_CLIMATOLOGY_OVERLOADS = .false.
  integer, public :: INPUT_VERSION_STRING = 0     ! Sub_rosa index
  integer, public :: OUTPUT_VERSION_STRING = 0    ! Sub_rosa index
  integer, public :: VERSION_COMMENT = 0          ! Sub_rosa index

!---------------------------- RCS Ident Info ---------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!-----------------------------------------------------------------------

contains

! subroutine SET_GLOBAL_SETTINGS ( ROOT, ForwardModelConfigDatabase ) !??? Restore when l2load isn't needed
  subroutine SET_GLOBAL_SETTINGS ( ROOT, ForwardModelConfigDatabase, FMC )
    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(ForwardModelConfig_T), dimension(:), pointer :: FORWARDMODELCONFIGDATABASE


!??? Begin temporary stuff to start up the forward model
  type(fwd_mdl_config) :: FMC
  integer :: IER
!??? End of temporary stuff to start up the forward model

    integer :: GSON                !??? Temporary for l2load
    character(len=255) :: LINE     !??? Temporary for l2load

    integer :: I         ! Index of son of root
    integer :: SON       ! Son of root

    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    do i = 2, nsons(root)-1 ! Skip names at beginning and end of section
      son = subtree(i,root)
      if ( node_id(son) == n_equal ) then
        select case ( decoration(subtree(1,son)) )
        case ( p_allow_climatology_overloads )
          allow_climatology_overloads = decoration(subtree(2,son)) == l_true
        case ( p_input_version_string )
          input_version_string = sub_rosa(subtree(2,son))
        case ( p_output_version_string )
          output_version_string = sub_rosa(subtree(2,son))
        case ( p_version_comment )
          version_comment = sub_rosa(subtree(2,son))
        end select
      else
        if ( node_id(son) == n_named ) son = subtree(2,son)
        select case ( get_spec_id(son) )
        case ( s_forwardModelGlobal ) !??? Begin temporary stuff for l2load
          call forwardModelGlobalSetup ( son )
        case ( s_forwardModel )
          call decorate (son, AddForwardModelConfigToDatabase ( &
            & forwardModelConfigDatabase, ConstructForwardModelConfig ( son ) ) )
        case ( s_l2load ) !??? More temporary stuff for l2load
          ! The only allowed field is the required ZVI field
          gson = subtree(2,son)
          call get_string ( sub_rosa(subtree(2,gson)), line ) ! ZVI file
          fmc%z = line(2:len_trim(line)-1)
          call l2_load ( fmc, ier=ier )
          !??? End temporary stuff for l2load
        end select
      end if
    end do

    if ( toggle(gen) ) call trace_end ( 'SET_GLOBAL_SETTINGS' )

  end subroutine SET_GLOBAL_SETTINGS

end module GLOBAL_SETTINGS

! $Log$
! Revision 2.8  2001/03/17 00:57:36  livesey
! Removed dump.
!
! Revision 2.7  2001/03/17 00:45:38  livesey
! Added ForwardModelConfigDatabase
!
! Revision 2.6  2001/03/09 02:30:13  vsnyder
! Allocate correct size for FMI and TFMI
!
! Revision 2.5  2001/03/09 00:24:30  vsnyder
! Do subscripts right
!
! Revision 2.4  2001/03/08 03:23:09  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.3  2001/03/07 22:46:05  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.2  2000/11/16 01:53:40  vsnyder
! Take timing back out.  Don't do it in sections that are only parameter settings.
!
! Revision 2.1  2000/11/16 01:45:25  vsnyder
! Implement timing.
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
