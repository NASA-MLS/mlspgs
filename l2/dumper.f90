! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DUMPER

! Dump various stuff so we can look at it.

  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: DUMP
  implicit NONE
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

! dump                   Dump the various user-defined types of mlsl2  
! === (end of toc) ===                                                   
! === (start of api) ===
! dump (type arg ) or     
! dump (type args(:), [int details] )      
!    where arg(s) may be among the following types:
!      { MLSChunk_t, hGrid_T, QuantityTemplates }      
! === (end of api) ===

  public :: DUMP

!---------------------------- RCS Ident Info ---------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!-----------------------------------------------------------------------

  interface DUMP
    module procedure DUMP_CHUNKS
    module procedure DUMP_a_HGRID
    module procedure DUMP_HGRIDS
  end interface

contains ! =====     Private Procedures     ============================
  ! ------------------------------------------------  DUMP_CHUNKS  -----
  subroutine DUMP_CHUNKS ( CHUNKS )

    use MLSCommon, only: MLSCHUNK_T

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

  ! ------------------------------------------------  DUMP_A_HGRID  -----
  subroutine DUMP_a_HGRID ( aHGRID )
    use HGridsDatabase, only: HGRID_T
    use STRING_TABLE, only: DISPLAY_STRING
    type(hGrid_T), intent(in) :: aHGRID
    integer :: J
      call output ( 'Name = ' )
      call display_string ( aHgrid%name )
      call output ( aHgrid%noProfs, before=' noProfs = ' )
      call output ( aHgrid%noProfsLowerOverlap, before=' lowerOverlap = ' )
      call output ( aHgrid%noProfsUpperOverlap, before=' upperOverlap = ', advance='yes' )
      call output ( ' prof       phi       geodLat           lon' )
      call output ( '          time     solarTime   solarZenith' )
      call output ( '      losAngle', advance='yes' )
      do j = 1, aHgrid%noProfs
        call output ( j, places=5 )
        call output ( aHgrid%phi(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%geodLat(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%lon(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%time(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarTime(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarZenith(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%losAngle(1,j), '(1x,1pg13.6)', advance='yes' )
      end do
  end subroutine DUMP_a_HGRID

  ! ------------------------------------------------  DUMP_HGRIDS  -----
  subroutine DUMP_HGRIDS ( HGRIDS )
    use HGridsDatabase, only: HGRID_T
    type(hGrid_T), intent(in) :: HGRIDS(:)
    integer :: I
    call output ( size(hgrids), before='HGRIDS: SIZE = ', advance='yes' )
    do i = 1, size(hgrids)
      call output ( i, 4, after=': ' )
      call dump ( hgrids(i) )
    end do
  end subroutine DUMP_HGRIDS

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DUMPER

! $Log$
! Revision 2.21  2004/05/18 01:07:33  vsnyder
! Repair broken Dump_a_HGrid and Dump_HGrids routines
!
! Revision 2.20  2004/05/01 04:03:47  vsnyder
! Get Dump_Quantity_Templates from QuantityTemplates instead of duplicating it
!
! Revision 2.19  2003/08/27 20:06:16  livesey
! More minor changes to dumping of chunks.
!
! Revision 2.18  2003/08/26 18:04:52  livesey
! Minor changes in dumping of chunks.
!
! Revision 2.17  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.16  2003/05/05 23:00:34  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.15.2.1  2003/03/27 23:17:10  vsnyder
! Use DUMP_a_HGRID in Dump_Hgrids instead of duplicating it
!
! Revision 2.15  2002/11/06 00:19:49  pwagner
! New chunk size info
!
! Revision 2.14  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.13  2002/08/22 01:22:20  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.12  2001/10/17 20:50:30  dwu
! a fix of cloud retrieval
!
! Revision 2.11  2001/08/06 18:37:36  pwagner
! Added Copyright statement
!
! Revision 2.10  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.9  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.8  2001/03/28 03:03:38  vsnyder
! Remove use, only's that aren't used
!
! Revision 2.7  2001/03/28 01:25:38  vsnyder
! Move DUMP_VGRIDS from dumper.f90 to VGrid.f90
!
