module GLOBAL_SETTINGS

  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT
  use MLSCommon, only: R8
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE

!??? Begin temporary stuff to start up the forward model
  use ForwardModelInterface, only: ForwardModelGlobalSetup, ForwardModelInfo_T
  use INIT_TABLES_MODULE, only: F_BILL, F_ZVI, S_ForwardModelGlobal, S_L2LOAD
  use L2_Load_M, only: L2_Load
  use L2_test_structures_m, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use String_Table, only: Get_String
  use TREE, only: NODE_ID
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

! subroutine SET_GLOBAL_SETTINGS ( ROOT, FwdModelInfo ) !??? Restore when l2load isn't needed
  subroutine SET_GLOBAL_SETTINGS ( ROOT, FwdModelInfo, FMC, FMI, TFMI )
    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(forwardModelInfo_T), intent(inout) :: FwdModelInfo ! From ForwardModelSetup

!??? Begin temporary stuff to start up the forward model
  type(fwd_mdl_config) :: FMC
  type(fwd_mdl_info), dimension(:), pointer :: FMI
  type(temporary_fwd_mdl_info), dimension(:), pointer :: TFMI
  integer :: IER
!??? End of temporary stuff to start up the forward model

    integer :: GSON                !??? Temporary for l2load
    integer :: J, K                !??? Temporary for l2load
    character(len=255) :: LINE     !??? Temporary for l2load

    integer :: I         ! Index of son of root
    integer :: SON       ! Son of root

    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    fwdModelInfo = forwardModelInfo_T(.false., .false., .false., 0.0_r8, &
      & .false., .false.)
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
          call forwardModelGlobalSetup ( son, fwdModelInfo )
          fmc%atmos_der = fwdModelInfo%atmos_der
          fmc%do_conv = fwdModelInfo%do_conv
          fmc%do_frqavg = fwdModelInfo%do_Freq_Avg
          fmc%spect_Der = fwdModelInfo%spect_Der
          fmc%temp_Der = fwdModelInfo%temp_Der
          fmc%zfrq = fwdModelInfo%the_Freq
        case ( s_l2load ) !??? More temporary stuff for l2load
          do j = 2, nsons(son)
            gson = subtree(j,son)
            if ( get_field_id(gson) == f_zvi ) then
              call get_string ( sub_rosa(subtree(2,gson)), line ) ! ZVI file
              fmc%z = line(2:len_trim(line)-1)
              call l2_load ( fmc, ier=ier )
          exit
            end if
          end do
          do j = 2, nsons(son)
            gson = subtree(j,son)
            if ( get_field_id(gson) == f_bill ) then
              if ( associated(fmi) ) deallocate ( fmi, stat=ier )
              if ( associated(tfmi) ) deallocate ( tfmi, stat=ier )
              allocate ( fmi(nsons(son)), stat=ier )
              if ( ier /= 0 ) call MLSMessage ( MLSmsg_Error, moduleName, &
                & MLSmsg_allocate // "fmi" )
              allocate ( tfmi(nsons(son)), stat=ier )
              if ( ier /= 0 ) call MLSMessage ( MLSmsg_Error, moduleName, &
                & MLSmsg_allocate // "tfmi" )
              do k = 2, nsons(gson)
                call get_string ( sub_rosa(subtree(k,gson)), line ) ! Bill file
                fmc%b = line(2:len_trim(line)-1)
                call l2_load ( fmc, fmi(k), tfmi(k), ier=ier )
              end do
            end if
          end do
          !??? End temporary stuff for l2load
        end select
      end if
    end do

    if ( toggle(gen) ) call trace_end ( 'SET_GLOBAL_SETTINGS' )

  end subroutine SET_GLOBAL_SETTINGS

end module GLOBAL_SETTINGS

! $Log$
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
