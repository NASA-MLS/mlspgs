module DUMPER

! Dump various stuff so we can look at it.

  use HGRID, only: HGRID_T
  use INIT_TABLES_MODULE, only: LIT_INDICES, PHYQ_INDICES
  use MLSCommon, only: MLSCHUNK_T
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use STRING_TABLE, only: DISPLAY_STRING
  use VGRID, only: VGRID_T
  implicit NONE
  private

  public :: DUMP

!---------------------------- RCS Ident Info ---------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!-----------------------------------------------------------------------

  interface DUMP
    module procedure DUMP_CHUNKS
    module procedure DUMP_HGRIDS
    module procedure DUMP_QUANTITY_TEMPLATES
    module procedure DUMP_VGRIDS
    module procedure DUMP_1D_DOUBLE
    module procedure DUMP_1D_INTEGER
    module procedure DUMP_2D_DOUBLE
    module procedure DUMP_2D_INTEGER
    module procedure DUMP_3D_DOUBLE
  end interface

  character, parameter :: AfterSub = '#'
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
      if ( .not. quantity_templates(i)%firstIndexChannel ) call output ( ' not' )
      call output ( ' firstIndexChannel' )
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
      if ( associated(quantity_templates(i)%signal) ) then
        ! ??? Dump the signal
      end if
      if ( quantity_templates(i)%radiometerIndex + &
        &  quantity_templates(i)%molecule /= 0 ) &
        &  call output ( '     ' )
      if ( quantity_templates(i)%radiometerIndex /= 0 ) then
        call output ( ' Radiometer = ' )
        call display_string ( quantity_templates(i)%radiometerIndex )
      end if
      if ( quantity_templates(i)%molecule /= 0 ) then
        call output ( ' Molecule = ' )
        call display_string ( lit_indices(quantity_templates(i)%molecule) )
      end if
      if ( quantity_templates(i)%radiometerIndex + &
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

  ! ------------------------------------------------  DUMP_VGRIDS  -----
  subroutine DUMP_VGRIDS ( VGRIDS )
    type(vGrid_T), intent(in) :: VGRIDS(:)
    integer :: I
    call output ( 'VGRIDS: SIZE = ' )
    call output ( size(vgrids), advance='yes' )
    do i = 1, size(vgrids)
      call output ( i, 4 )
      call output ( ': Name = ' )
      call display_string ( vgrids(i)%name )
      call output ( ' noSurfs = ' )
      call output ( vgrids(i)%noSurfs )
      call output ( ' verticalCoordinate = ' )
      call display_string ( lit_indices(vgrids(i)%verticalCoordinate) )
      call dump ( vgrids(i)%surfs, ' Surfs = ' )
    end do
  end subroutine DUMP_VGRIDS

  ! ---------------------------------------------  DUMP_1D_DOUBLE  -----
  subroutine DUMP_1D_DOUBLE ( ARRAY, NAME )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer :: J, K
    if ( size(array) == 1 ) then
      if ( present(name) ) call output ( name )
      call output ( array(1), '(1x,1pg13.6)', advance='yes' )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do j = 1, size(array), 5
        call output ( j, 4 ); call output ( afterSub )
        do k = j, min(j+4, size(array))
          call output ( array(k), '(1x,1pg13.6)' )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer :: J, K
    if ( size(array) == 1 ) then
      if ( present(name) ) call output ( name )
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do j = 1, size(array), 10
        call output ( j, 4 ); call output ( afterSub )
        do k = j, min(j+9, size(array))
          call output ( array(k), 6 )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_INTEGER

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
  subroutine DUMP_2D_DOUBLE ( ARRAY, NAME )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer :: I, J, K
    if ( size(array) == 1 ) then
      if ( present(name) ) call output ( name )
      call output ( array(1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 5
          call output ( i, 4 )
          call output ( j, 4 ); call output ( afterSub )
          do k = j, min(j+4, size(array,2))
            call output ( array(i,k), '(1x,1pg13.6)' )
          end do
          call output ( '', advance='yes' )
        end do
      end do
    end if
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  subroutine DUMP_2D_INTEGER ( ARRAY, NAME )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer :: I, J, K
    if ( size(array) == 1 ) then
      if ( present(name) ) call output ( name )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          call output ( i, 4 )
          call output ( j, 4 ); call output ( afterSub )
          do k = j, min(j+9, size(array,2))
            call output ( array(i,k), 6 )
          end do
          call output ( '', advance='yes' )
        end do
      end do
    end if
  end subroutine DUMP_2D_INTEGER

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    integer :: I, J, K, L
    if ( size(array) == 1 ) then
      if ( present(name) ) call output ( name )
      call output ( array(1,1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            call output ( i, 4 )
            call output ( j, 4 )
            call output ( k, 4 ); call output ( afterSub )
            do l = k, min(k+4, size(array,3))
              call output ( array(i,j,l), '(1x,1pg13.6)' )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end do
    end if
  end subroutine DUMP_3D_DOUBLE
end module DUMPER

! $Log,v $
