! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
!=============================================================================

! Set up the forward model.  Interface from the retrieve step to the
! forward model.

!??? Do we want a forward model database ???

  use Init_Tables_Module, only: field_first, field_indices, field_last
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decorate, Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, &
    & Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: Vector_T

  implicit NONE
  private
  public :: ForwardModel, ForwardModelInfo_T, ForwardModelSetup

  type ForwardModelInfo_T
    integer :: Foo !??? Just because the compiler insists
  end type ForwardModelInfo_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

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

    integer :: Error                    ! Error level seen so far
    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    integer :: Son                      ! Some subtree of root.

    ! Error message codes
    integer, parameter :: Twice = 1     ! A field appears twice

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
      field = decoration(subtree(1,son))
      if ( got(field) ) call announceError ( twice, field )
      got(field) = .true.
      select case ( field )
      case default
        ! Shouldn't get here if the type checker worked
      end select
    end do ! i = 2, nsons(key)
    if ( toggle(gen) ) call trace_end ( "ForwardModelSetup" )

    contains
    ! --------------------------------------------  AnnounceError  -----
    subroutine AnnounceError ( Code, FieldIndex )
      integer, intent(in) :: Code       ! Index of error message
      integer, intent(in) :: FieldIndex ! f_...

      integer :: Source

      error = max(error,1)
      source = source_ref ( son )
      call output ( 'At line '  )
      call output ( mod(source,256) )
      call output ( ', column ' )
      call output ( source/256 )
      call output ( ' ForwardModelSetup complained: ' )
      select case ( code )
      case ( twice )
        call output ( 'the field ' )
        call display_string ( field_indices(fieldIndex) )
        call output ( ' shall not appear twice.', advance='yes' )
      end select
    end subroutine AnnounceError
  end subroutine ForwardModelSetup

  ! -----------------------------------------------  ForwardModel  -----
  subroutine ForwardModel ( FwdModelInfo, FwdModelIn, Jacobian, F, RowBlock, &
    &                       FwdModelOut )
    type(forwardModelInfo_T), intent(in) :: FwdModelInfo ! From ForwardModelSetup
    type(vector_T), intent(in) :: FwdModelIn           ! ???
    type(matrix_T), intent(inout), optional :: Jacobian
    type(vector_T), intent(in), optional :: F          ! Computed radiances
    integer, intent(in), optional :: RowBlock          ! With which block of
    ! rows of F and Jacobian are we computing? All of them if absent.
    type(vector_T), intent(inout), optional :: FwdModelOut  ! ???
  end subroutine ForwardModel

end module ForwardModelInterface

! $Log$
! Revision 2.2  2001/02/08 00:56:11  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.1  2001/02/07 00:52:27  vsnyder
! Initial commit
!
