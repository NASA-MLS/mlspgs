! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DestroyCommand_m

! Destroy vectors and matrices in the databases

  implicit None
  private

  public :: DestroyCommand

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DestroyCommand ( Key, Matrices, Vectors )

    use Init_Tables_Module, only: F_AllMatrices, F_AllVectors, F_Matrix, F_Vector
    use MatrixModule_1, only: DestroyMatrix, Dump, Matrix_Database_T
    use MoreTree, only: Get_Boolean, Get_Field_ID
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Toggles, only: Toggle, Gen
    use Trace_m, only: Trace_begin, Trace_end
    use Tree, only: Decoration, NSons, Subtree
    use VectorsModule, only: DestroyVectorInfo, Dump, Vector_T

    integer, intent(in) :: Key ! Root of parse subtree
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (vector_T), dimension(:), pointer :: Vectors

    logical :: DEEBUG = .false.
    integer :: J, K, MatrixToKill, Son, SourceVectorIndex

    if ( toggle(gen) ) call trace_begin ( 'DestroyCommand' )
    ! Here we're to try to shrink the vector database by destroying a vector
    ! or the matrix database by destroying a matrix
    ! Loop over the instructions
    do j = 2, nsons(key)
      son = subtree(j,key)  ! The argument
      select case ( get_field_id(son) )
      case ( f_allMatrices )
        if ( get_boolean(son) ) then
          do k = 1, size(matrices)
            call DestroyMatrix ( matrices(k) )
          end do
        end if
      case ( f_allVectors )
        if ( get_boolean(son) ) then
          do k = 1, size(vectors)
            call DestroyVectorInfo ( vectors(k) )
          end do
        end if
      case ( f_matrix )
        do k = 2, nsons(son)
          matrixToKill = decoration(decoration(subtree(k,son)))
          if ( DEEBUG ) then
            ! if ( matrices(matrixToKill)%matrix%name /= 0 ) then
             ! call output ( '   Matrix Name = ' )
             ! call display_string ( matrices(matrixToKill)%matrix%name )
             call dump ( matrices(matrixToKill), -1 )
            ! end if
          end if

          call DestroyMatrix ( matrices(matrixToKill) )
        end do
      case ( f_vector )
        do k = 2, nsons(son)
          sourceVectorIndex = decoration(decoration(subtree(k,son)))
          if ( DEEBUG ) then
            if ( vectors(sourceVectorIndex)%name /= 0 ) then
              call display_string ( vectors(sourceVectorIndex)%name, &
                & before='   Vector Name = ' )
            else
              call output ( '  ' )
            end if
            if ( vectors(sourceVectorIndex)%template%name /= 0 ) then
              call display_string ( vectors(sourceVectorIndex)%template%name, &
                & before=' Template_Name = ', advance='yes' )
            end if
            call output ( ' -- vector database before removal --', advance='yes' )
            call dump ( vectors, details=-2 )
          end if

          call DestroyVectorInfo ( vectors(sourceVectorIndex) )
     !    vectorindex = rmVectorFromDatabase ( vectors, vectors(sourceVectorIndex) )
          if ( DEEBUG ) then
            call output ( ' -- vector database after removal --', advance='yes' )
            call dump ( vectors, details=-2 )
          end if
        end do
      case default ! Can't get here if type checker worked
      end select
    end do
    if ( toggle(gen) ) call trace_end ( 'DestroyCommand' )

  end subroutine DestroyCommand

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here


end module DestroyCommand_m

! $Log$
! Revision 2.2  2006/08/03 20:06:56  vsnyder
! Added /allvectors and /allmatrices
!
! Revision 2.1  2006/08/02 19:51:43  vsnyder
! Move from Fill module
!
