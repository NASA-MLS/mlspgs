module DUMPER

! Dump various stuff so we can look at it.

  use DUMP_0, only: AfterSub, DUMP
  use HGRID, only: HGRID_T
  use INIT_TABLES_MODULE, only: LIT_INDICES, PHYQ_INDICES
  use MLSCommon, only: MLSCHUNK_T
  use MLSSignals_m, only: radiometers, signals, DUMP, GetRadiometerName, GetModuleName
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use STRING_TABLE, only: DISPLAY_STRING
  use VGRID, only: VGRID_T
  implicit NONE
  private

  public :: DUMP

!---------------------------- RCS Ident Info ---------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
!-----------------------------------------------------------------------

  interface DUMP
    module procedure DUMP_CHUNKS
    module procedure DUMP_HGRIDS
    module procedure DUMP_QUANTITY_TEMPLATES
  end interface

contains ! =====     Private Procedures     ============================
  ! ------------------------------------------------  DUMP_CHUNKS  -----
  subroutine DUMP_CHUNKS ( CHUNKS )
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
      call output ( '      accumulatedMAFs = ' )
      call output ( chunks(i)%accumulatedMAFs, advance='yes' )
    end do
  end subroutine DUMP_CHUNKS

  ! ------------------------------------------------  DUMP_HGRIDS  -----
  subroutine DUMP_HGRIDS ( HGRIDS )
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
        call output ( hgrids(i)%phi(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%geodLat(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%lon(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%time(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%solarTime(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%solarZenith(j), '(1x,1pg13.6)' )
        call output ( hgrids(i)%losAngle(j), '(1x,1pg13.6)', advance='yes' )
      end do
    end do
  end subroutine DUMP_HGRIDS

  ! ------------------------------------  DUMP_QUANTITY_TEMPLATES  -----
  subroutine DUMP_QUANTITY_TEMPLATES ( QUANTITY_TEMPLATES )
    type(QuantityTemplate_T), intent(in) :: QUANTITY_TEMPLATES(:)
    integer :: I
    character (len=80) :: str
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
      call output ( quantity_templates(i)%InstanceLen )
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
      if ( quantity_templates(i)%radiometer + &
        &  quantity_templates(i)%molecule /= 0 ) &
        &  call output ( '     ' )
      if ( quantity_templates(i)%radiometer /= 0 ) then
        call output ( ' Radiometer = ' )
        call GetRadiometerName ( quantity_templates(i)%radiometer, str )
        call output ( str )
      end if
      if ( quantity_templates(i)%instrumentModule /= 0 ) then
        call output ( ' Instrument Module = ' )
        call GetModuleName ( quantity_templates(i)%instrumentModule, str )
        call output ( str )
      end if
      if ( quantity_templates(i)%molecule /= 0 ) then
        call output ( ' Molecule = ' )
        call display_string ( lit_indices(quantity_templates(i)%molecule) )
      end if
      call output ( '', advance = 'yes')
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
    end do
  end subroutine DUMP_QUANTITY_TEMPLATES
end module DUMPER

! $Log$
