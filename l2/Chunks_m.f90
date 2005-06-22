! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Chunks_m
!=============================================================================

  implicit NONE
  private
  public :: Dump, MLSChunk_T

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  type MLSChunk_T
    integer :: firstMAFIndex   ! Index of first MAF in the chunk
    integer :: lastMAFIndex    ! Index of last MAF in the chunk
    integer :: noMAFsLowerOverlap ! Number of MAFs in the lower overlap region
    integer :: noMAFsUpperOverlap ! Number of MAFs in the upper overlap region
    integer :: chunkNumber              ! Index of this chunk
    integer, dimension(:), pointer :: HGridOffsets => NULL()
    ! This for each chunk is the index of the first non-overlapped profile in 
    ! each hGrid into the relevant output (l2gp?) file.
    integer, dimension(:), pointer :: HGridTotals => NULL()
    ! This is somewhat repetetive.  It's the total number of profiles in
    ! the output hGrid.  It's only really used in parallel runs.
  end type MLSChunk_T

  interface DUMP
    module procedure Dump_Chunks
  end interface

contains ! =====     Private Procedures     ============================

  ! ------------------------------------------------  DUMP_CHUNKS  -----
  subroutine DUMP_CHUNKS ( CHUNKS )

    use Output_M, only: Output

    type(MLSChunk_t), intent(in) :: CHUNKS(:)
    integer :: I
    call output ( 'CHUNKS: SIZE = ' )
    call output ( size(chunks), advance='yes' )
    do i = 1, size(chunks)
      call output ( i, before=' Chunk ', advance='yes' )
      call output ( chunks(i)%firstMAFIndex, before='  firstMAFIndex: ' )
      call output ( chunks(i)%lastMAFIndex, before='  lastMAFIndex: ', advance='yes' )
      call output ( chunks(i)%noMAFsLowerOverlap, before='  noMAFsLowerOverlap: ' )
      call output ( chunks(i)%noMAFsUpperOverlap, before='  noMAFsUpperOverlap: ', advance='yes' )
      call output ( chunks(i)%firstMAFIndex + chunks(i)%noMAFsLowerOverlap, &
        & before='  1st non-overlapped MAF: ' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%noMAFsUpperOverlap, &
        & before='  last non-overlapped MAF: ', advance='yes' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%firstMAFIndex + 1, &
        & before='  chunk size: ' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%firstMAFIndex &
        & - chunks(i)%noMAFsUpperOverlap - chunks(i)%noMAFsLowerOverlap + 1, &
        & before='  non-overlapped chunk size: ', advance='yes' )
    end do
  end subroutine DUMP_CHUNKS

!=======================================================================
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Chunks_m

! $Log$
! Revision 2.2  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2004/05/19 19:15:07  vsnyder
! Initial commit
!
