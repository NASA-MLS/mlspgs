! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Species_Data_M

  use ForwardModelVectorTools, only: QtyStuff_T
  use MLSCommon, only: RP
  use SpectroscopyCatalog_m, only: Catalog_t

  implicit NONE
  private
  public :: Get_Species_Data, Destroy_Beta_Group, Destroy_Catalog_Extract
  public :: Destroy_Species_Data, Dump

! Beta group type declaration:
  type, public :: Beta_Group_T
    integer :: N_Elements
    integer, pointer  :: Cat_Index(:) => NULL()
    real(rp), pointer :: Ratio(:) => NULL()
  end type Beta_Group_T

! Species stuff
  type, public :: SPS_T
    integer :: No_Mol ! Size of several arrays == # molecule groups + # PFA
    integer :: No_LBL ! Number of Line-by-line molecule groups, <= no_mol
    type(beta_group_t), pointer :: Beta_Group(:) => NULL() ! 1:no_lbl
    type(catalog_t), pointer :: Catalog(:,:) => NULL()     ! sidebands,1:no_species
    integer, pointer :: Mol_Cat_Index(:) => NULL()         ! 1:no_mol; indices
      ! of elements of FwdModelConf%molecules that are positive.
    type(qtyStuff_t), pointer :: Qtys(:) => NULL()         ! 1:no_mol
  end type SPS_T

  interface Dump
    module procedure Dump_Beta_Group, Dump_Sps_Data
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!-----------------------------------------------------------------------
contains

  ! -------------------------------------------  Get_Species_Data  -----
  subroutine Get_Species_Data ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    & Radiometer, Sps )

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: LIT_INDICES, L_ISOTOPERATIO, L_VMR
    use MLSCommon, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
      & MLSMSG_Warning
    use MLSSets, only: FINDFIRST
    use MLSSignals_m, only: GetRadiometerFromSignal
    use SpectroscopyCatalog_m, only: Catalog, Dump, Empty_Cat, &
      & Line_t, Lines
    use String_table, only: GET_STRING
    use Toggles, only: Switches
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T

    ! Inputs

    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    integer, intent(in) :: Radiometer

    ! Outputs

    type (sps_t), intent(out) :: Sps

    ! Local variables

    real(rp) :: Beta_Ratio
    logical :: DoThis                   ! Flag for lines in catalog item
    type (VectorValue_T), pointer :: F  ! An arbitrary species
    integer :: I, IER, J, K, L
    integer, dimension(:), pointer :: LINEFLAG ! /= 0 => Use this line
    ! (noLines per species)
    integer, dimension(size(fwdModelConf%molecules,1)+1) :: Molecules_Temp
    character (len=32) :: molName       ! Name of a molecule
    integer :: NLines                   ! count(lineFlag)
    integer :: NoMol                    ! size(fwdModelConf%molecules) -- includes
      ! negative ones that indicate molecule grouping
    integer  :: NoNonPFA                ! No. of non-PFA Molecules under
      ! consideration --  includes negative ones that indicate molecule grouping
    integer :: Polarized                ! -1 => One of the selected lines is Zeeman split
    ! +1 => None of the selected lines is Zeeman split
    integer :: SIGIND                   ! Signal index, loop counter
    integer :: SV_I
    integer :: S                        ! Sideband index
    type (catalog_T), pointer :: thisCatalogEntry
    type (line_T), pointer :: thisLine

    if ( .not. associated ( catalog ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'No spectroscopy catalog has been defined' )

    nullify ( lineFlag )

    noMol = size(fwdModelConf%molecules)
    noNonPFA = fwdModelConf%firstPFA - 1
    sps%no_lbl = count(fwdModelConf%molecules(1:noNonPFA) > 0)
    sps%no_mol = sps%no_lbl + noMol - fwdModelConf%firstPFA + 1

    if ( sps%no_mol == 0 ) return

    allocate ( sps%beta_group(sps%no_lbl), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'sps%beta_group' )

    call allocate_test ( sps%mol_cat_index, sps%no_mol, 'sps%mol_cat_index', moduleName )

    k = max(1,noNonPFA-sps%no_lbl)
    do i = 1, sps%no_lbl
      call allocate_test ( sps%beta_group(i)%cat_index, k, 'beta_group%cat_index', moduleName )
      call allocate_test ( sps%beta_group(i)%ratio, k, 'beta_group%ratio', moduleName )
      sps%beta_group(i)%n_elements = 0
      sps%beta_group(i)%ratio = 0.0
      sps%beta_group(i)%cat_index = 0
    end do

    if ( noNonPFA == sps%no_lbl ) then ! No grouping, .eqv. sps%no_mol == noMol

      ! Work out beta group etc for line-by-line
      sv_i = 0
      beta_ratio = 1.0_rp   ! Always, for single element (no grouping)
      do j = 1, sps%no_lbl
        l = fwdModelConf%molecules(j)
        !        if ( l == l_extinction ) CYCLE
        sv_i = sv_i + 1
        sps%mol_cat_index(sv_i) = j
        sps%beta_group(sv_i)%n_elements   = 1
        sps%beta_group(sv_i)%cat_index(1) = j
        sps%beta_group(sv_i)%ratio(1)     = beta_ratio
      end do

      ! All that's needed for PFA is sps%mol_cat_index(sv_i), to get indices
      ! for Beta_Path and dBetaD*
      do j = 1, noNonPFA+1, noMol
        l = fwdModelConf%molecules(j)
        !        if ( l == l_extinction ) CYCLE
        sv_i = sv_i + 1
        sps%mol_cat_index(sv_i) = j
      end do

    else

      molecules_temp(1:noMol) = fwdModelConf%molecules(1:noMol)
      molecules_temp(noMol+1) = noMol ! Anything positive will do

      sps%mol_cat_index = PACK((/(i,i=1,noMol)/), fwdModelConf%molecules > 0)

      ! Work out beta group etc for line-by-line
      sv_i = 0
      do j = 1, noNonPFA
        k = molecules_temp(j)
        l = abs(k)
        !        if ( l == l_extinction ) CYCLE
        beta_ratio = 1.0_rp
        if ( k > 0 ) then
          if ( molecules_temp(j+1) > 0 ) then
            sv_i = sv_i + 1
            sps%beta_group(sv_i)%n_elements   = 1
            sps%beta_group(sv_i)%cat_index(1) = j
            sps%beta_group(sv_i)%ratio(1)     = beta_ratio
          end if
        else
          if ( molecules_temp(j-1) > 0) sv_i = sv_i + 1
          f => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_isotopeRatio, molecule=l, noError=.TRUE., &
            & config=fwdModelConf )
          if ( associated ( f ) ) beta_ratio = f%values(1,1)
          i = sps%beta_group(sv_i)%n_elements + 1
          sps%beta_group(sv_i)%n_elements   = i
          sps%beta_group(sv_i)%cat_index(i) = j
          sps%beta_group(sv_i)%ratio(i)     = beta_ratio
        end if
      end do

    end if

    ! Work out which spectroscopy we're going to need ------------------------

    allocate ( &
      & sps%catalog(fwdModelConf%sidebandStart:fwdModelConf%sidebandStop,noNonPFA), &
      & stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'sps%catalog' )

    sps%catalog = empty_cat
    do s = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
      j = 0
      do sv_i = 1, noNonPFA
        j = j + 1
        ! Skip if the next molecule is negative (indicates that this one is a
        ! parent)
        if ( (sv_i < noNonPFA) .and. (fwdModelConf%molecules(sv_i) > 0) ) then
          if (fwdModelConf%molecules(sv_i+1) < 0 ) then
            ! sps%catalog springs into existence with %lines and %polarized null
            call allocate_test ( sps%catalog(s,j)%lines, 0, &
              & 'sps%catalog(?,?)%lines(0)', moduleName )
            call allocate_test ( sps%catalog(s,j)%polarized, 0, &
              & 'sps%catalog(?,?)%polarized(0)', moduleName )
            cycle
          end if
        end if
        l = abs(fwdModelConf%molecules(sv_i))
        thisCatalogEntry => Catalog(FindFirst(catalog%molecule, l ) )
        sps%catalog(s,j) = thisCatalogEntry
        ! Don't deallocate them by mistake -- sps%catalog is a shallow copy
        nullify ( sps%catalog(s,j)%lines, sps%catalog(s,j)%polarized )
        if ( associated ( thisCatalogEntry%lines ) ) then
          ! Now subset the lines according to the signal we're using
          call allocate_test ( lineFlag, size(thisCatalogEntry%lines), &
            &  'lineFlag', moduleName )
          lineFlag = 0
          if ( fwdModelConf%allLinesInCatalog ) then
            ! NOTE: If allLinesInCatalog is set, then no lines can be polarized,
            ! this is checked for in ForwardModelSupport.
            lineFlag = 1
          else
            do k = 1, size ( thisCatalogEntry%lines )
              thisLine => lines(thisCatalogEntry%lines(k))
              if ( associated(thisLine%signals) ) then
                polarized = 1 ! not polarized
                ! Work out whether to do this line
                do sigInd = 1, size(fwdModelConf%signals)
                  if ( fwdModelConf%allLinesForRadiometer ) then
                    doThis = .false.
                    do i = 1, size(thisLine%signals)
                      ! Tried to make GetRadiometerFromSignal elemental, but compile time
                      ! in LF95 (optimized) for Construct.f90 went nuts! :-(
                      if ( GetRadiometerFromSignal ( thisLine%signals(i) ) == &
                        & fwdModelConf%signals(sigInd)%radiometer ) then
                        doThis = .true.
                        if ( .not. fwdModelConf%polarized ) &
                          exit   ! loop over signals for line -- no need to check for
                        ! polarized lines
                        if ( associated(thisLine%polarized) ) then
                          if ( thisLine%polarized(i) ) then
                            polarized = -1 ! polarized
                            exit   ! loop over signals for line -- one signal that sees a
                            ! polarized line is enough to turn on the polarized
                            ! method
                          end if
                        end if
                      end if
                    end do ! End loop over signals for line
                  else
                    ! Not doing all lines for radiometer, be more selective
                    doThis = any ( &
                      & ( thisLine%signals == fwdModelConf%signals(sigInd)%index ) .and. &
                      & ( ( thisLine%sidebands == 0 ) .or. ( thisLine%sidebands == s ) ) )
                    if ( fwdModelConf%polarized .and. doThis .and. &
                      & associated(thisLine%polarized) ) then
                      if ( any(thisLine%polarized) ) polarized = -1 ! polarized
                    end if
                  end if
                  
                  if ( fwdModelConf%sidebandStart == fwdModelConf%sidebandStop ) &
                    & doThis = doThis .and. &
                    & any( ( thisLine%sidebands == fwdModelConf%sidebandStart ) .or. &
                    & ( thisLine%sidebands == 0 ) )
                  if ( doThis ) then
                    lineFlag(k) = polarized
                    if ( polarized < 0 .or. .not. fwdModelConf%polarized ) &
                      exit   ! loop over signals requested in fwm
                  end if
                end do ! End loop over signals requested in fwm
              end if
            end do     ! End loop over lines
          end if       ! End case where allLinesInCatalog not set

          ! Check we have at least one line for this specie

          nLines = count(lineFlag /= 0)
          if ( nLines == 0 .and. all ( sps%catalog(s,j)%continuum == 0.0 ) &
            & .and. (index(switches, '0sl') > 0) ) then
            call get_string ( lit_indices(l), molName )
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'No relevant lines or continuum for '//trim(molName) )
          end if
          call allocate_test ( sps%catalog(s,j)%lines, nLines, &
            & 'sps%catalog(?,?)%lines', moduleName )
          call allocate_test ( sps%catalog(s,j)%polarized, nLines, &
            & 'sps%catalog(?,?)%polarized', moduleName )
          sps%catalog(s,j)%lines = pack ( thisCatalogEntry%lines, lineFlag /= 0 )
          sps%catalog(s,j)%polarized = pack ( lineFlag < 0, lineFlag /= 0 )
          call deallocate_test ( lineFlag, 'lineFlag', moduleName )

        else

          ! No lines for this species.  However, its continuum is still valid 
          ! so don't set it to empty.
          ! Won't bother checking that continuum /= 0 as if it was then
          ! presumably having no continuum and no lines it wouldn't be in the catalog!
          call allocate_test ( sps%catalog(s,j)%lines, 0, 'sps%catalog(?,?)%lines(0)', &
            & moduleName )
          call allocate_test ( sps%catalog(s,j)%polarized, 0, 'sps%catalog(?,?)%polarized(0)', &
            & moduleName )
        end if
      end do         ! Loop over species
    end do                              ! Loop over sidebands

    ! Get state vector quantities for species
    allocate ( sps%qtys(sps%no_mol), stat = i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'sps%qtys' )

    do sv_i = 1 , sps%no_mol
      sps%qtys(sv_i)%qty => GetQuantityForForwardModel(fwdmodelin, fwdmodelextra, &
        &  quantitytype=l_vmr, molIndex=sps%mol_cat_index(sv_i), config=fwdModelConf, &
        &  radiometer=radiometer, foundInFirst=sps%qtys(sv_i)%foundInFirst )
    end do

    if ( index(switches,'sps') /= 0 ) call dump ( sps )

  end subroutine Get_Species_Data

  ! -----------------------------------------  Destroy_Beta_Group  -----
  subroutine Destroy_Beta_Group ( Beta_Group )

  ! Destroy the catalog extract prepared by Get_Species_Data

    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type(beta_group_t), pointer :: Beta_Group(:)

    integer :: I

    do i = 1, size(beta_group)
      call deallocate_test ( beta_group(i)%cat_index, 'beta_group(i)%cat_index', &
        & moduleName )
      call deallocate_test ( beta_group(i)%ratio, 'beta_group(i)%ratio', &
        & moduleName )
    end do

    deallocate ( beta_group, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'beta_group' )

  end subroutine Destroy_Beta_Group

  ! ------------------------------------  Destroy_Catalog_Extract  -----
  subroutine Destroy_Catalog_Extract ( My_Catalog )

  ! Destroy the catalog extract prepared by Get_Species_Data

    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type(catalog_t), pointer :: My_Catalog(:,:)

    integer :: I, J
    do j = lbound(my_catalog,1), ubound(my_catalog,1), 2
      do i = 1, size(my_catalog,2)
        ! Note that we don't deallocate the signals/sidebands stuff for each line
        ! as these are shallow copies of the main spectroscopy catalog stuff
        call deallocate_test ( my_catalog(j,i)%lines, 'my_catalog(?,?)%lines', &
          & moduleName )
        call deallocate_test ( my_catalog(j,i)%polarized, 'my_catalog(?,?)%polarized', &
          & moduleName )
      end do
    end do

    deallocate ( my_catalog, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'My_Catalog' )

  end subroutine Destroy_Catalog_Extract

  ! ---------------------------------------  Destroy_Species_Data  -----
  subroutine Destroy_Species_Data ( Species_Data )

  ! Destroy the SPS_T structure

    use Allocate_Deallocate, only: DEALLOCATE_TEST
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type(sps_t), intent(inout) :: Species_Data

    integer :: I

    if ( species_data%no_mol == 0 ) return
    call destroy_beta_group ( species_data%beta_group )
    call destroy_catalog_extract ( species_data%catalog )
    call deallocate_test ( species_data%mol_cat_index, 'Mol_Cat_Index', moduleName )
    deallocate ( species_data%qtys, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'Species_Data%Qtys' )

  end subroutine Destroy_Species_Data

  ! --------------------------------------------  Dump_Beta_Group  -----
  subroutine Dump_Beta_Group ( Beta_Group, Name )

    use Dump_0, only: Dump
    use Output_m, only: Output

    type(beta_group_t), intent(in) :: Beta_Group(:)
    character(len=*), intent(in), optional :: Name

    integer :: I

    call output ( 'Beta group' )
    if ( present(name) ) call output ( ' ' // trim(name) )
    call output ( ', SIZE = ' )
    call output ( size(beta_group), advance='yes' )
    do i = 1, size(beta_group)
      call output ( 'Item ' )
      call output ( i, advance='yes' )
      call dump ( beta_group(i)%cat_index(:beta_group(i)%n_elements), name='Cat_Index' )
      call dump ( beta_group(i)%ratio(:beta_group(i)%n_elements), name='Ratio' )
    end do
  end subroutine Dump_Beta_Group

  ! ----------------------------------------------  Dump_Sps_Data  -----
  subroutine Dump_Sps_Data ( Sps, Details )

    use Dump_0, only: Dump
    use Output_m, only: NewLine, Output
    use SpectroscopyCatalog_m, only: Dump
    use String_Table, only: Display_String

    type(sps_t), intent(in) :: Sps
    integer, intent(in), optional :: Details ! for dumping spectroscopy catalog

    integer :: I

    call output ( 'Species data', advance='yes' )
    call output ( sps%no_mol, before=' Total molecules = ' )
    call output ( sps%no_lbl, before=', Line-by-line molecules = ', advance='yes' )
    call dump ( sps%beta_group )
    call dump ( sps%catalog, details=details )
    call dump ( sps%mol_cat_index(:sps%no_lbl), name='Mol_Cat_Index (non-PFA)' )
    call dump ( sps%mol_cat_index(sps%no_lbl+1:sps%no_mol), &
      & name='Mol_Cat_Index (PFA)', lbound=sps%no_lbl+1 )
    call output ( 'QtyStuff:', advance='yes' )
    do i = 1, size(sps%qtys)
      call output ( ' Quantity ' )
      if ( sps%qtys(i)%qty%template%name /= 0 ) then
        call display_string ( sps%qtys(i)%qty%template%name )
      else
        call output ( '???' )
      end if
      if ( sps%qtys(i)%foundInFirst ) call output ( ', Found in first' )
      call newLine
    end do

  end subroutine Dump_Sps_Data

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module  Get_Species_Data_M

! $Log$
! Revision 2.12  2004/08/03 22:06:45  vsnyder
! Inching further toward PFA
!
! Revision 2.11  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.10  2004/06/10 00:59:56  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.9  2004/03/22 18:23:56  livesey
! Added handling of AllLinesInCatalog flag (precludes polarized)
!
! Revision 2.8  2003/10/09 23:32:31  pwagner
! SIPS version should stop complaining about nolines
!
! Revision 2.7  2003/07/15 18:17:04  livesey
! Catalog now split by sideband
!
! Revision 2.6  2003/06/27 00:59:53  vsnyder
! Simplify interface to Get_Species_Data
!
! Revision 2.5  2003/05/24 02:25:53  vsnyder
! Set the polarized flag correctly -- well, at least differently
!
! Revision 2.4  2003/05/21 22:15:36  vsnyder
! Dump my_catalog and beta_group if the 'bgrp' switch is set
!
! Revision 2.3  2003/05/17 01:19:32  vsnyder
! Remove unreferenced USE name, futzing
!
! Revision 2.2  2003/05/16 23:48:59  livesey
! Removed references to spectags (note old code had a bug when looking for
! h2o_r??, which led to a 0.08K error in band 2).
!
! Revision 2.1  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 1.1.2.9  2003/05/01 23:53:03  livesey
! Bug fix was being overzelous with setting my_catalog(j)=empty_cat
!
! Revision 1.1.2.8  2003/03/22 04:03:45  vsnyder
! Move Beta_Group_T and Dump_Beta_Group from get_beta_path to Get_Species_Data.
! Write Destroy_Beta_Group.
!
! Revision 1.1.2.7  2003/03/13 02:03:09  vsnyder
! Initialize some uninitialized variables
!
! Revision 1.1.2.6  2003/03/05 03:27:08  vsnyder
! Don't clobber spectroscopy catalog by way of shallow copy
!
! Revision 1.1.2.5  2003/03/01 03:18:39  vsnyder
! Fix bugs in calculation of the 'polarized' field
!
! Revision 1.1.2.4  2003/02/27 23:21:47  vsnyder
! Put polarized flag in my_catalog.  Add Destroy_Species_Data subroutine.
!
! Revision 1.1.2.3  2003/02/21 21:04:30  vsnyder
! Just to make CVS happy about a merge that didn't do anything
!
! Revision 1.1.2.2  2003/02/18 22:58:37  pwagner
! Compatible with FullForwardModel
!
