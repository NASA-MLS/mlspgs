! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DUMPER

! Dump various stuff so we can look at it.

  use OUTPUT_M, only: OUTPUT
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
    module procedure DUMP_QUANTITY_TEMPLATES
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
      call output ( i, 4 )
      call output ( ': firstMAFIndex = ' )
      call output ( chunks(i)%firstMAFIndex )
      call output ( ' lastMAFIndex = ' )
      call output ( chunks(i)%lastMAFIndex, advance='yes' )
      call output ( '      noMAFsLowerOverlap = ' )
      call output ( chunks(i)%noMAFsLowerOverlap )
      call output ( ' noMAFsUpperOverlap = ' )
      call output ( chunks(i)%noMAFsUpperOverlap, advance='yes' )
      call output ( '      1st non-overlap chunk = ' )
      call output ( chunks(i)%firstMAFIndex + chunks(i)%noMAFsLowerOverlap )
      call output ( '      last non-overlap chunk = ' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%noMAFsUpperOverlap, &
        & advance='yes' )
      call output ( '      chunk size= ' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%firstMAFIndex )
      call output ( '      non-overlap chunk size= ' )
      call output ( chunks(i)%lastMAFIndex - chunks(i)%firstMAFIndex &
        & - chunks(i)%noMAFsUpperOverlap - chunks(i)%noMAFsLowerOverlap + 1)
      call output ( '      accumulatedMAFs = ' )
      call output ( chunks(i)%accumulatedMAFs, advance='yes' )
    end do
  end subroutine DUMP_CHUNKS

  ! ------------------------------------------------  DUMP_A_HGRID  -----
  subroutine DUMP_a_HGRID ( aHGRID )
    use HGRID, only: HGRID_T
    type(hGrid_T), intent(in) :: aHGRID
    integer :: J
      do j = 1, ahgrid%noProfs
        call output ( ahgrid%phi(j), '(1x,1pg13.6)' )
        call output ( ahgrid%geodLat(j), '(1x,1pg13.6)' )
        call output ( ahgrid%lon(j), '(1x,1pg13.6)' )
        call output ( ahgrid%time(j), '(1x,1pg13.6)' )
        call output ( ahgrid%solarTime(j), '(1x,1pg13.6)' )
        call output ( ahgrid%solarZenith(j), '(1x,1pg13.6)' )
        call output ( ahgrid%losAngle(j), '(1x,1pg13.6)', advance='yes' )
      end do
  end subroutine DUMP_a_HGRID

  ! ------------------------------------------------  DUMP_HGRIDS  -----
  subroutine DUMP_HGRIDS ( HGRIDS )
    use HGRID, only: HGRID_T
    use STRING_TABLE, only: DISPLAY_STRING
    type(hGrid_T), intent(in) :: HGRIDS(:)
    integer :: I, J
    call output ( 'HGRIDS: SIZE = ' )
    call output ( size(hgrids), advance='yes' )
    do i = 1, size(hgrids)
      call output ( i, 4 )
      call output ( ': Name = ' )
      call display_string ( hgrids(i)%name )
      call output ( ' noProfs = ' )
      call output ( hgrids(i)%noProfs, advance='yes' )
      call output ( '      noProfsLowerOverlap = ' )
      call output ( hgrids(i)%noProfsLowerOverlap )
      call output ( ' noProfsUpperOverlap = ' )
      call output ( hgrids(i)%noProfsUpperOverlap, advance='yes' )
      call output ( ' prof          phi       geodLat           lon' )
      call output ( '          time     solarTime   solarZenith' )
      call output ( '      losAngle', advance='yes' )
      do j = 1, hgrids(i)%noProfs
        call output ( j, 4 )
        call dump ( hgrids(i) )
      end do
    end do
  end subroutine DUMP_HGRIDS

  ! ------------------------------------  DUMP_QUANTITY_TEMPLATES  -----
  subroutine DUMP_QUANTITY_TEMPLATES ( QUANTITY_TEMPLATES, DETAILS )

    use DUMP_0, only: DUMP
    use Intrinsic, only: LIT_INDICES, PHYQ_INDICES
    use MLSSignals_m, only: signals, DUMP, GetRadiometerName, GetModuleName
    use QuantityTemplates, only: QuantityTemplate_T
    use STRING_TABLE, only: DISPLAY_STRING

    type(QuantityTemplate_T), intent(in) :: QUANTITY_TEMPLATES(:)
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
    !                                        ! >0   => Do dump arrays
    !                                        ! Default 1
    integer :: I, MyDetails
    character (len=80) :: Str
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ( 'QUANTITY_TEMPLATES: SIZE = ' )
    call output ( size(quantity_templates), advance='yes' )
    do i = 1, size(quantity_templates)
      call output ( i, 4 )
      call output ( ': Name = ' )
      call display_string ( quantity_templates(i)%name )
      call output ( ' Id = ' ); call output ( quantity_templates(i)%id )
      call output ( ' quantityType = ' )
      call display_string ( lit_indices(quantity_templates(i)%quantityType), &
        & advance='yes' )
      call output ( '      NoInstances = ' )
      call output ( quantity_templates(i)%noInstances )
      call output ( ' NoSurfs = ' )
      call output ( quantity_templates(i)%noSurfs )
      call output ( ' noChans = ' )
      call output ( quantity_templates(i)%noChans, advance='yes' )
      call output ( '      ' )
      if ( .not. quantity_templates(i)%coherent ) call output ( 'in' )
      call output ( 'coherent ' )
      if ( .not. quantity_templates(i)%stacked ) call output ( 'non' )
      call output ( 'stacked ' )
      if ( .not. quantity_templates(i)%regular ) call output ( 'ir' )
      call output ( 'regular ' )
      if ( quantity_templates(i)%logBasis ) then
        call output ('log-')
      else
        call output ('linear-')
      endif
      call output ('basis ' )  
      if ( .not. quantity_templates(i)%minorFrame ) call output ( 'non' )
      call output ( 'minorFrame', advance='yes' )
      call output ( '      NoInstancesLowerOverlap = ' )
      call output ( quantity_templates(i)%noInstancesLowerOverlap )
      call output ( ' NoInstancesUpperOverlap = ' )
      call output ( quantity_templates(i)%noInstancesUpperOverlap, advance='yes' )
      call output ( '      BadValue = ' )
      call output ( quantity_templates(i)%badValue )
      call output ( ' Unit = ' )
      call display_string ( phyq_indices(quantity_templates(i)%unit) )
      call output ( ' ScaleFactor = ' )
      call output ( quantity_templates(i)%scaleFactor, advance='yes' )
      call output ( '      InstanceLen = ' )
      call output ( quantity_templates(i)%InstanceLen, advance='yes' )
      if ( myDetails < 0 ) then
        call dump ( quantity_templates(i)%surfs, '  Surfs = ' )
        call dump ( quantity_templates(i)%phi, '      Phi = ' )
        call dump ( quantity_templates(i)%geodLat, '      GeodLat = ' )
        call dump ( quantity_templates(i)%lon, '      Lon = ' )
        call dump ( quantity_templates(i)%time, '      Time = ' )
        call dump ( quantity_templates(i)%solarTime, '      SolarTime = ' )
        call dump ( quantity_templates(i)%solarZenith, '      SolarZenith = ' )
        call dump ( quantity_templates(i)%losAngle, '      LosAngle = ' )
        if ( associated(quantity_templates(i)%mafIndex) ) then
          call dump ( quantity_templates(i)%mafIndex, '      MAFIndex = ' )
          call dump ( quantity_templates(i)%mafCounter, '      MAFCounter = ' )
        end if
        if ( associated(quantity_templates(i)%frequencies) ) then
          call output ( '      FrequencyCoordinate = ' )
          call output ( quantity_templates(i)%frequencyCoordinate )
          call dump ( quantity_templates(i)%frequencies, ' Frequencies = ' )
        end if
      end if
      if ( quantity_templates(i)%radiometer /= 0 ) then
        call output ( '      Radiometer = ' )
        call GetRadiometerName ( quantity_templates(i)%radiometer, str )
        call output ( trim(str), advance='yes' )
      end if
      if ( quantity_templates(i)%molecule + &
        &  quantity_templates(i)%instrumentModule /= 0 ) then
        call output ( '     ' )
        if ( quantity_templates(i)%molecule /= 0 ) then
          call output ( ' Molecule = ' )
          call display_string ( lit_indices(quantity_templates(i)%molecule) )
        end if
        if ( quantity_templates(i)%instrumentModule /= 0 ) then
          call output ( ' Instrument Module = ' )
          call GetModuleName ( quantity_templates(i)%instrumentModule, str )
          call output ( trim(str) )
        end if
        call output ( '', advance = 'yes')
      end if
      if ( myDetails > 0 ) then
        if ( quantity_templates(i)%signal /= 0 ) then
          call dump ( signals( (/ quantity_templates(i)%signal /) ) )
        end if
        if ( quantity_templates(i)%radiometer + &
          &  quantity_templates(i)%molecule /= 0 ) &
          &  call output ( '', advance='yes' )
        if ( associated(quantity_templates(i)%surfIndex) ) then
          call dump ( quantity_templates(i)%surfIndex, '      SurfIndex = ' )
        end if
        if ( associated(quantity_templates(i)%chanIndex) ) then
          call dump ( quantity_templates(i)%chanIndex, '      ChanIndex = ' )
        end if
      end if
    end do
  end subroutine DUMP_QUANTITY_TEMPLATES
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DUMPER

! $Log$
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
