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

  public :: DESTROYCOMMAND

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DestroyCommand ( Key, Matrices, Vectors, Grids )

    use GRIDDEDDATA, only: GRIDDEDDATA_T, &
      & DESTROYGRIDDEDDATA, DESTROYGRIDDEDDATADATABASE
    use INIT_TABLES_MODULE, only: F_ALLGRIDDEDDATA, F_ALLMATRICES, &
      & F_ALLVECTORS, F_BOOLEAN, F_GRID, F_MATRIX, &
      & F_VECTOR
    use MATRIXMODULE_1, only: DESTROYMATRIXDATABASE, DESTROYMATRIX, DUMP, &
      & MATRIX_DATABASE_T
    use MLSL2OPTIONS, only: REMOVERUNTIMEBOOLEAN
    use MLSSTRINGS, only: LOWERCASE
    use MORETREE, only: GET_BOOLEAN, GET_FIELD_ID
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING, GET_STRING
    use TOGGLES, only: TOGGLE, GEN
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE
    use VECTORSMODULE, only: DESTROYVECTORDATABASE, DESTROYVECTORINFO, DUMP, &
      & VECTOR_T

    integer, intent(in) :: Key ! Root of parse subtree
    type (matrix_database_T), dimension(:), pointer :: Matrices
    type (vector_T), dimension(:), pointer :: Vectors
    type (griddedData_T), dimension(:), pointer :: Grids
    ! local variables
    character(len=80) :: BOOLEANNAME    ! E.g., 'BQTYS'
    logical :: DEEBUG = .false.
    integer :: gridID, J, K, MatrixToKill, Son, SourceVectorIndex
    integer :: gson
    integer :: Me = -1                  ! String index for trace

    ! Executable

    call trace_begin ( me, 'DestroyCommand', key, cond=toggle(gen) )
    ! Here we're to try to shrink the vector database by destroying a vector
    ! or the matrix database by destroying a matrix
    ! Loop over the instructions
    booleanName = ' '
    do j = 2, nsons(key)
      son = subtree(j,key)  ! The argument
      gson = son
      if ( nsons(gson) > 1 ) gson = subtree(2,gson)
      select case ( get_field_id(son) )
      case ( f_allGriddedData )
        if ( get_boolean(son) ) call DestroyGriddedDataDatabase ( grids )
      case ( f_allMatrices )
        if ( get_boolean(son) ) call destroyMatrixDatabase ( matrices )
      case ( f_allVectors )
        if ( get_boolean(son) ) call destroyVectorDatabase ( vectors )
      case (f_Boolean)
        call get_string ( sub_rosa(gson), booleanName, strip=.true. )
        booleanName = lowerCase( booleanName )
        call removeRunTimeBoolean( booleanName )
      case ( f_grid )
        do k = 2, nsons(son)
          gridID = decoration(decoration(subtree(k,son)))
          call DestroyGriddedData ( grids(gridID) )
        end do
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
    call trace_end ( 'DestroyCommand', cond=toggle(gen) )

  end subroutine DestroyCommand

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------


end module DestroyCommand_m

! $Log$
! Revision 2.7  2013/08/30 02:45:36  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.6  2013/05/22 20:22:26  pwagner
! May destroy r/t Booleans
!
! Revision 2.5  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2008/09/19 23:54:31  pwagner
! May now Destroy GriddedData
!
! Revision 2.3  2006/08/04 18:08:23  vsnyder
! Simplify /allMatrices and /allVectors
!
! Revision 2.2  2006/08/03 20:06:56  vsnyder
! Added /allvectors and /allmatrices
!
! Revision 2.1  2006/08/02 19:51:43  vsnyder
! Move from Fill module
!
