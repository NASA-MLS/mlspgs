module GLOBAL_SETTINGS

  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT, &
    & S_TIME
  use OUTPUT_M, only: OUTPUT
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE

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

  subroutine SET_GLOBAL_SETTINGS ( ROOT )
    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree

    integer :: I         ! Index of son of root
    integer :: SON       ! Son of root
    double precision :: T1, T2     ! for timing
    logical :: TIMING

    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    timing = .false.     ! timing isn't on.
    do i = 2, nsons(root)-1 ! Skip names at beginning and end of section
      son = subtree(i,root)
      select case ( decoration(subtree(1,son)) )
      case ( p_allow_climatology_overloads )
        allow_climatology_overloads = decoration(subtree(2,son)) == l_true
      case ( p_input_version_string )
        input_version_string = sub_rosa(subtree(2,son))
      case ( p_output_version_string )
        output_version_string = sub_rosa(subtree(2,son))
      case ( p_version_comment )
        version_comment = sub_rosa(subtree(2,son))
      case ( s_time )
        if ( timing ) then
          call cpu_time ( t2 )
          call output ( "Timing for Set_Global_Settings = " )
          call output ( t2-t1, advance = 'yes' )
          timing = .false.
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select
    end do
    if ( timing ) then
      call cpu_time ( t2 )
      call output ( "Timing for Set_Global_Settings = " )
      call output ( t2-t1, advance = 'yes' )
    end if

  if ( toggle(gen) ) call trace_end ( 'SET_GLOBAL_SETTINGS' )

  end subroutine SET_GLOBAL_SETTINGS

end module GLOBAL_SETTINGS

! $Log$
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
