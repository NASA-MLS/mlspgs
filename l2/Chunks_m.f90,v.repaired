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

  use Dates_Module, only: TAI93S2HID
  use MLSCommon, only: MLSChunk_T
  implicit none
  private
  public :: Dump, Initialize, MLSChunk_T

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface Dump
    module procedure Dump_Chunks, Dump_One_Chunk
  end interface

contains ! =====     Private Procedures     ============================

  ! ------------------------------------------------  Dump_Chunks  -----
  subroutine Dump_Chunks ( Chunk )
    ! Let's write them as a table instead of Dumping each one indivdually

    use HighOutput, only: Tab, OutputNamedValue
    use Output_m, only: Blanks, NewLine, Output

    type(MLSChunk_t), intent(in) :: Chunk(:)
    integer :: I
    character(len=*), parameter :: The80 = &
      & '12345678901234567890123456789012345678901234567890123456789012345678901234567890'
    character(len=*), parameter :: TheDecades = &
      & '         1         2         3         4         5         6         7         8'
    call outputNamedValue ( 'Size(chunks)', size(chunk), options='--Headline' )
    ! call output( The80, advance='yes' )
    ! call output( TheDecades, advance='yes' )
    ! Column headers: 3 lines
    !                       MAF 
    !            Index      Overlaps  Non-overlap        chunk non-overlap      phi
    ! Chunk   First  Last  Lower Upper    First Last      size    size     start   end
    call tab ( 3 )
    call output( 'MAF', advance='yes' )
    call tab ( 3 )
    call output( 'Index', advance='no' )
    call tab ( 7 )
    call output( 'Overlaps', advance='no' )
    call tab ( 10 )
    call output( 'Non-overlaps', advance='no' )
    call tab ( 14 )
    call output( 'chunk', advance='no' )
    call blanks ( 1 )
    call output( 'Non-overlap', advance='no' )
    call tab ( 18 )
    call output( 'phi', advance='no' )
    call newLine
    
    call output( 'chunk', advance='no' )
    call tab
    call output( 'First', advance='no' )
    call tab
    call output( 'Last', advance='no' )
    call tab
    call tab
    call output( 'Lower', advance='no' )
    call tab
    call output( 'Upper', advance='no' )
    call tab
    call output( 'First', advance='no' )
    call tab
    call output( 'Last', advance='no' )
    call tab
    call tab
    call output( 'size', advance='no' )
    call tab
    call output( 'size', advance='no' )
    call tab
    call tab
    call output( 'start', advance='no' )
    call tab
    call output( 'end', advance='no' )
    call newLine
    do i = 1, size(chunk)
      call output ( i, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%firstMAFIndex, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%lastMAFIndex, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%noMAFsLowerOverlap, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%noMAFsUpperOverlap, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%firstMAFIndex + chunk(i)%noMAFsLowerOverlap, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%lastMAFIndex - chunk(i)%noMAFsUpperOverlap, advance='no' )
      call tab
      call tab
      call output ( chunk(i)%lastMAFIndex - chunk(i)%firstMAFIndex+1, advance='no' )
      call tab
      call output ( chunk(i)%lastMAFIndex - chunk(i)%firstMAFIndex &
      & - chunk(i)%noMAFsUpperOverlap - chunk(i)%noMAFsLowerOverlap + 1, advance='no' )
      call tab
      call output ( chunk(i)%phiStart, advance='no', format='(F10.2)' )
      call blanks( 1 )
      call output ( chunk(i)%phiEnd, advance='no', format='(F10.2)' )
      call tab
      call output ( tai93s2hid(chunk(i)%StartTime, leapsec=.true.), advance='no', format='(F10.2)' )
      call blanks( 1 )
      call output ( tai93s2hid(chunk(i)%EndTime, leapsec=.true.), advance='no', format='(F10.2)' )
      call NewLine
    end do
  end subroutine Dump_Chunks

  subroutine Dump_Chunks_old ( Chunks )
    ! Man, this one was terrible

    use Output_m, only: Output

    type(MLSChunk_t), intent(in) :: Chunks(:)
    integer :: I
    call output ( 'CHUNKS: SIZE = ' )
    call output ( size(chunks), advance='yes' )
    do i = 1, size(chunks)
      call output ( i, before=' Chunk ', advance='yes' )
      call dump(chunks(i))
    end do
  end subroutine Dump_Chunks_old

  ! ---------------------------------------------  Dump_One_Chunk  -----
  subroutine Dump_One_Chunk ( Chunk )

    use Dump_0, only: Dump
    use Output_m, only: Output

    type(MLSChunk_t), intent(in) :: Chunk
    if ( chunk%chunkNumber > -1 ) call output ( chunk%chunkNumber, before='  chunkNumber: ' )
    call output ( chunk%firstMAFIndex, before='  firstMAFIndex: ' )
    call output ( chunk%lastMAFIndex, before='  lastMAFIndex: ', advance='yes' )
    call output ( chunk%noMAFsLowerOverlap, before='  noMAFsLowerOverlap: ' )
    call output ( chunk%noMAFsUpperOverlap, before='  noMAFsUpperOverlap: ', &
      & advance='yes' )
    call output ( chunk%firstMAFIndex + chunk%noMAFsLowerOverlap, &
      & before='  1st non-overlapped MAF: ' )
    call output ( chunk%lastMAFIndex - chunk%noMAFsUpperOverlap, &
      & before='  last non-overlapped MAF: ', advance='yes' )
    call output ( chunk%lastMAFIndex - chunk%firstMAFIndex + 1, &
      & before='  chunk size: ' )
    call output ( chunk%lastMAFIndex - chunk%firstMAFIndex &
      & - chunk%noMAFsUpperOverlap - chunk%noMAFsLowerOverlap + 1, &
      & before='  non-overlapped chunk size: ', advance='yes' )
    call output ( '  phi start, end: ' )
    call output ( (/chunk%phiStart, chunk%phiEnd /), format='(F10.2)', &
      & advance='no' )
    call output ( '  hid start, end: ' )
    call output ( (/tai93s2hid(chunk%StartTime, leapsec=.true.), &
      & tai93s2hid(chunk%EndTime, leapsec=.true.) /), format='(F10.2)', &
      & advance='yes' )
    if ( associated(chunk%HGridOffsets) ) &
      & call dump( chunk%HGridOffsets, 'HGrid offsets' )
    if ( associated(chunk%HGridTotals) ) &
      & call dump( chunk%HGridTotals, 'HGrid Totals' )
  end subroutine Dump_One_Chunk

  ! ------------------------------------------------  Initialize  -----
  ! Override each chunk component's default values if given a precursor
  subroutine Initialize ( Chunk, Precursor )
    type(MLSChunk_T)                       :: Chunk
    type(MLSChunk_T), optional, intent(in) :: Precursor
    if ( .not. present(precursor) ) return
    chunk%abandoned            = precursor%abandoned
    chunk%firstMAFIndex        = precursor%lastMAFIndex + 1
    chunk%lastMAFIndex         = chunk%firstMAFIndex + precursor%lastMAFIndex - precursor%firstMAFIndex
    chunk%noMAFsLowerOverlap   = precursor%noMAFsLowerOverlap
    chunk%noMAFsUpperOverlap   = precursor%noMAFsUpperOverlap
    chunk%chunkNumber          = precursor%chunkNumber       
    chunk%phiStart             = precursor%phiEnd      
    chunk%phiEnd               = 2*precursor%phiEnd - precursor%phiStart
    chunk%StartTime            = precursor%EndTime      
    chunk%EndTime              = 2*precursor%EndTime - precursor%StartTime
  end subroutine Initialize

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Chunks_m

! $Log$
! Revision 2.13  2020/07/09 23:52:01  pwagner
! Added Start,EndTime components to MLSChunk_T
!
! Revision 2.12  2018/11/01 23:17:33  pwagner
! Dont need to print column numbers any more
!
! Revision 2.11  2018/03/02 00:59:44  pwagner
! Improve how ChunkDivide appears when Dump-ed
!
! Revision 2.10  2016/07/28 00:41:41  vsnyder
! Remove unreferenced USE
!
! Revision 2.9  2015/06/19 00:36:39  pwagner
! Moved MLSChunk_T to MLSCommon in lib
!
! Revision 2.8  2012/06/21 00:40:28  pwagner
! Added phi start and end to be used someday by HGrid
!
! Revision 2.7  2011/06/29 21:53:55  pwagner
! Added initialize, improved dump
!
! Revision 2.6  2010/05/23 03:34:29  honghanh
! Initialize firstMafIndex, lastMAfIndex and a couple of other variables
!
! Revision 2.5  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2009/02/03 17:58:38  pwagner
! Added the abandoned field
!
! Revision 2.3  2005/09/21 23:22:42  pwagner
! Dump may dump a single chunk
!
! Revision 2.2  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2004/05/19 19:15:07  vsnyder
! Initial commit
!
