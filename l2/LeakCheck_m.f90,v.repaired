! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module LeakCheck_m

  ! Scan the l2cf looking for vector, matrix and destroy commands.
  ! Report all vectors and matrices that are created but not destroyed.

  implicit NONE
  private
  public :: LeakCheck

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine LeakCheck ( Root )
    use init_tables_module, only: f_matrix, f_vector, &
      & s_destroy, s_matrix, s_vector
    use MoreTree, only: Get_Field_ID, Get_Spec_ID, StartErrorMessage
    use Output_m, only: Output
    use String_Table, only: Display_String, String_Table_Size
    use Tree, only: Node_id, Nsons, Sub_Rosa, Subtree, Where_At => Where
    use Tree_Types, only: N_CF, N_Named, N_Spec_Args
    integer, intent(in) :: Root ! of the parse tree

    type Record_t
      character(len=6) :: What  ! matrix, vector
      integer :: Name  ! Sub_Rosa
      integer :: Define=0, Destroy=0 ! Tree index, for Where_At
    end type Record_t

    integer :: GSon, GGSon, I, J, K, Key, L, Name, NumRec, Son

    type(record_t) :: Record(string_Table_Size()) ! Certainly big enough

    numRec = 0
    do i = 1, nsons(root)
      son = subtree(i, root)
      if ( node_id(son) /= n_cf ) cycle ! Skip type etc. definitions
      do j = 2, nsons(son) - 1      ! skip names at begin/end of section
        gson = subtree(j, son)
        if ( node_id(gson) == n_named ) then
          key = subtree(2, gson)
          name = sub_rosa(subtree(1,gson))
        else
          key = gson
          name = 0
        end if
        if ( node_id(key) /= n_spec_args ) cycle ! Skip parameter settings
        select case ( get_spec_id(key) )
        case ( s_destroy )
          do k = 2, nsons(key)
            ggson = subtree(k, key)
            select case ( get_field_id(ggson) )
            case ( f_matrix )
              do l = 2, nsons(ggson)
                call delete ( 'matrix', sub_rosa(subtree(l,ggson)), ggson )
              end do
            case ( f_vector )
              do l = 2, nsons(ggson)
                call delete ( 'vector', sub_rosa(subtree(l,ggson)), ggson )
              end do
            end select
          end do
        case ( s_matrix )
          call add ( 'matrix', name, key )
        case ( s_vector )
          call add ( 'vector', name, key )
        end select
      end do
    end do

    do i = 1, numRec
      if ( record(i)%destroy == 0 ) then ! Found a leak
        call startErrorMessage ( record(i)%define )
        call display_string ( record(i)%name, before=record(i)%what // ' ' )
        call output ( ' created but not destroyed.', advance='yes' )
      end if
    end do

  contains
    subroutine Add ( What, Who, Where )
      character(len=*), intent(in) :: What ! Vector, Matrix?
      integer, intent(in) :: Who, Where ! Sub_Rosa, root
      integer :: I
      do i = 1, numRec
        if ( who == record(i)%name ) then
          call startErrorMessage ( where )
          call display_string ( who, before='Duplicate ' // trim(what) // ' ' )
          call output ( ' defined.', advance='yes' )
        end if
      end do
      numRec = numRec + 1
      record(numRec) = record_t(what,who,where,0)
    end subroutine Add

    subroutine Delete ( What, Who, Where )
      use Lexer_core, only: Print_source
      use Output_m, only: NewLine
      character(len=*), intent(in) :: What ! Vector, Matrix?
      integer, intent(in) :: Who, Where ! Sub_Rosa, root
      integer :: I
      do i = 1, numRec
        if ( who == record(i)%name ) then
          if ( what /= record(i)%what ) then
            call startErrorMessage ( where )
            call display_string ( who )
            call output ( ' created at ' )
            call print_source ( where_at(record(i)%define) )
            call output (' as ' // record(i)%what // &
              &          ' but destroyed as ' // what, advance='yes' )
          end if
          if ( record(i)%destroy /= 0 ) then
            call startErrorMessage ( where )
            call display_string ( who, before=what // ' ' )
            call output ( ' already destroyed at ' )
            call print_source ( where_at(where) )
            call newLine
          end if
          record(i)%destroy = where
          return
        end if
      end do
      call startErrorMessage ( where )
      call display_string ( who )
      call output ( ' destroyed but never created.', advance='yes' )
    end subroutine Delete

  end subroutine LeakCheck

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LeakCheck_m

! $Log$
! Revision 2.4  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.3  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2006/07/27 03:48:11  vsnyder
! Skip parameter settings
!
! Revision 2.1  2006/07/27 03:09:28  vsnyder
! Initial commit
!
