! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Species_Data_m

  ! Get species data for the molecules in the beta groups.

  implicit NONE
  private
  public :: Get_Species_Data, Destroy_Species_Data

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)), save :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !-----------------------------------------------------------------------

contains

  subroutine Get_Species_Data ( FwdModelConf, FwdModelIn, FwdModelExtra )

    ! Fill in the Beta_Groups in FwdModelConf; create and fill fwdModelConf%catalog.
    ! DeriveFromForwardModelConfig needs to be called BEFORE this routine,
    ! because it allocates FwdModelConf%channels, which size is used here.

    use Allocate_Deallocate, only: Allocate_Test
    use ForwardModelConfig, only: Dump, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: LIT_INDICES, L_ISOTOPERATIO, L_NONE, L_VMR
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
      & MLSMSG_Warning
    use MLSSignals_m, only: GetRadiometerFromSignal
    use SpectroscopyCatalog_m, only: Catalog, Dump, Empty_Cat, Line_t, &
      & Lines, MostLines
    use String_table, only: GET_STRING
    use Toggles, only: Switches
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf ! Fills Beta_Group
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra

    integer :: B         ! Index for beta groups 
    integer :: C         ! Index for fwdModelConf%catalog, or size(fwdModelConf%channels)
    logical :: DoThis    ! Flag for lines in catalog item
    integer, save :: DumpFWM = -1
    type (VectorValue_T), pointer :: F  ! An arbitrary species
    integer :: I         ! Index for signals for a line
    integer :: K         ! Index in main spectroscopy catalog
    integer :: L         ! Index for lines, or number of lines
    integer, pointer :: LINEFLAG(:) ! /= 0 => Use this line
    integer :: M         ! Index for molecules in beta groups, or size thereof
    integer, target :: MaxLineFlag(mostLines)
    character(len=32) :: MolName
    integer :: N         ! Molecule name
    integer :: Polarized ! -1 => One of the selected lines is Zeeman split
                         ! +1 => None of the selected lines is Zeeman split
    integer :: S         ! Index for sidebands                
    integer :: STAT      ! Status from allocate or deallocate 
    type (line_T), pointer :: ThisLine
    integer :: Z         ! Index for fwdModelConf%Signals

    if ( dumpFWM < 0 ) then ! done only once
      if ( index(switches,'fwmg') > 0 ) dumpFWM = 1 ! Dump but don't stop
      if ( index(switches,'fwmG') > 0 ) dumpFWM = 2 ! Dump and stop
    end if

    ! Get isotope ratios for molecules in a beta group, else 1.0 if not a group
    do b = 1, size(fwdModelConf%beta_group)
      if ( fwdModelConf%beta_group(b)%group ) then ! A molecule group
        ! First LBL molecules' ratios
        do m = 1, size(fwdModelConf%beta_group(b)%lbl_molecules)
          f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_isotopeRatio, &
            & molecule=fwdModelConf%beta_group(b)%lbl_molecules(m), &
            & noError=.TRUE., config=fwdModelConf )
          fwdModelConf%beta_group(b)%lbl_ratio(m)     = 1.0
          if ( associated ( f ) ) & ! Have an isotope ratio
            & fwdModelConf%beta_group(b)%lbl_ratio(m) = f%values(1,1)
        end do ! m
        ! Now PFA molecules' ratios
        do m = 1, size(fwdModelConf%beta_group(b)%pfa_molecules)
          f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_isotopeRatio, &
            & molecule=fwdModelConf%beta_group(b)%pfa_molecules(m), &
            & noError=.TRUE., config=fwdModelConf )
          fwdModelConf%beta_group(b)%pfa_ratio(m)     = 1.0
          if ( associated ( f ) ) & ! Have an isotope ratio
            & fwdModelConf%beta_group(b)%pfa_ratio(m) = f%values(1,1)
        end do ! m
!     else ! Not a molecule group, but this is handled in ForwardModelSupport
!       fwdModelConf%beta_group(b)%lbl_ratio(1)     = 1.0
!       fwdModelConf%beta_group(b)%pfa_ratio(1)     = 1.0
      end if
    end do ! b

    ! Allocate the spectroscopy catalog extract

    c = 0
    do b = 1, size(fwdModelConf%beta_group) ! Get total catalog size
      c = c + size(fwdModelConf%beta_group(b)%lbl_molecules)
    end do

    allocate ( fwdModelConf%catalog(fwdModelConf%sidebandStart:fwdModelConf%sidebandStop,c), &
      & stat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Allocate//'fwdModelConf%catalog' )

    fwdModelConf%catalog = empty_cat

    ! Work out the spectroscopy we're going to need.
    do s = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2
      c = 0
      do b = 1, size(fwdModelConf%beta_group)
        do m = 1, size(fwdModelConf%beta_group(b)%lbl_molecules)
          c = c + 1
          fwdModelConf%beta_group(b)%cat_index(m) = c
          n = fwdModelConf%beta_group(b)%lbl_molecules(m)
          if ( catalog(n)%molecule == l_none ) then
            call get_string ( lit_indices(n), molName )
            call MLSMessage ( MLSMSG_Error, moduleName, &
              & 'No spectroscopy catalog for ' // molName )
          end if
          fwdModelConf%catalog(s,c) = catalog(n)
          ! Don't deallocate them by mistake -- fwdModelConf%catalog is a shallow copy
          nullify ( fwdModelConf%catalog(s,c)%lines, fwdModelConf%catalog(s,c)%polarized )
          if ( associated ( catalog(n)%lines ) ) then
            ! Now subset the lines according to the signal we're using
            lineFlag => MaxLineFlag(:size(catalog(n)%lines))
            lineFlag = 0
            if ( fwdModelConf%allLinesInCatalog ) then
              ! NOTE: If allLinesInCatalog is set, then no lines can be polarized;
              ! this is checked for in ForwardModelSupport.
              lineFlag = 1
            else
              do l = 1, size ( catalog(n)%lines )
                thisLine => lines(catalog(n)%lines(l))
                if ( associated(thisLine%signals) ) then
                  polarized = 1 ! not polarized
                  ! Work out whether to do this line
                  do z = 1, size(fwdModelConf%signals)
                    if ( fwdModelConf%allLinesForRadiometer ) then
                      doThis = .false.
                      do i = 1, size(thisLine%signals)
                        ! Tried to make GetRadiometerFromSignal elemental, but compile time
                        ! in LF95 (optimized) for Construct.f90 went nuts! :-(
                        if ( GetRadiometerFromSignal ( thisLine%signals(i) ) == &
                          & fwdModelConf%signals(z)%radiometer ) then
                          doThis = .true.
                          if ( .not. fwdModelConf%polarized ) &
                            exit   ! loop over signals for line -- no need to check for
                          ! polarized lines
                          if ( associated(thisLine%polarized) ) then
                            if ( thisLine%polarized(i) ) then
                              polarized = -1 ! polarized
                              exit   ! loop over signals for line -- one signal
                              ! that sees a polarized line is enough to turn on
                              ! the polarized method
                            end if
                          end if
                        end if
                      end do ! End loop over signals for line
                    else
                      ! Not doing all lines for radiometer, be more selective
                      doThis = any ( &
                        & ( thisLine%signals == fwdModelConf%signals(z)%index ) .and. &
                        & ( ( thisLine%sidebands == 0 ) .or. ( thisLine%sidebands == s ) ) )
                      if ( fwdModelConf%polarized .and. doThis .and. &
                        & associated(thisLine%polarized) ) then
                        if ( any(thisLine%polarized) ) polarized = -1 ! polarized
                      end if
                    end if

                    if ( fwdModelConf%sidebandStart == fwdModelConf%sidebandStop ) &
                      & doThis = doThis .and. &
                      & any( ( thisLine%sidebands == fwdModelConf%sidebandStart ) &
                      & .or. ( thisLine%sidebands == 0 ) )
                    if ( doThis ) then
                      lineFlag(l) = polarized
                      if ( polarized < 0 .or. .not. fwdModelConf%polarized ) &
                        exit   ! loop over signals requested in fwm
                    end if
                  end do ! z End loop over signals requested in fwm
                end if
              end do     ! l End loop over lines
            end if       ! End case where allLinesInCatalog not set

            ! Check we have at least one line for this specie

            l = count(lineFlag /= 0)
            if ( l == 0 .and. all ( fwdModelConf%catalog(s,c)%continuum == 0.0 ) &
              & .and. (index(switches, '0sl') > 0) ) then
              call get_string ( lit_indices(n), molName )
              call MLSMessage ( MLSMSG_Warning, ModuleName, &
                & 'No relevant lines or continuum for '//trim(molName) )
            end if
            call allocate_test ( fwdModelConf%catalog(s,c)%lines, l, &
              & 'fwdModelConf%catalog(?,?)%lines', moduleName )
            call allocate_test ( fwdModelConf%catalog(s,c)%polarized, l, &
              & 'fwdModelConf%catalog(?,?)%polarized', moduleName )
            fwdModelConf%catalog(s,c)%lines = pack ( catalog(n)%lines, lineFlag /= 0 )
            fwdModelConf%catalog(s,c)%polarized = pack ( lineFlag < 0, lineFlag /= 0 )

          else

            ! No lines for this species.  However, its continuum is still valid 
            ! so don't set it to empty.
            ! Won't bother checking that continuum /= 0 as if it were then
            ! presumably having no continuum and no lines it wouldn't be in the
            ! catalog!
            call allocate_test ( fwdModelConf%catalog(s,c)%lines, 0, &
              & 'fwdModelConf%catalog(?,?)%lines(0)', moduleName )
            call allocate_test ( fwdModelConf%catalog(s,c)%polarized, 0, &
              & 'fwdModelConf%catalog(?,?)%polarized(0)', moduleName )
          end if
        end do ! m Molecules in fwdModelConf
      end do ! b Beta groups
    end do ! s Sidebands

    ! Get state vector quantities for species
    do b = 1, size(fwdModelConf%beta_group)
      fwdModelConf%beta_group(b)%qty%qty => getQuantityForForwardModel ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_vmr, molIndex=b,    &
        &  config=fwdModelConf, radiometer=fwdModelConf%signals(1)%radiometer, &
        &  foundInFirst=fwdModelConf%beta_group(b)%qty%foundInFirst )
    end do ! b

    if ( dumpFWM > 0 ) then
      call dump ( fwdModelConf, moduleName )
      call dump ( fwdModelConf%catalog )
      if ( dumpFWM > 1 ) stop
    end if

  end subroutine Get_Species_Data

  ! -------------------------------------------  Destroy_Species_Data  -----
  subroutine Destroy_Species_Data ( FwdModelConf )

  ! Destroy the spectroscopy catalog extract that was allocated by
  ! Get_Species_Data.

    use Allocate_Deallocate, only: Deallocate_Test
    use ForwardModelConfig, only: ForwardModelConfig_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf

    integer :: S, C ! Sideband index or status, catalog index

    do s = lbound(fwdModelConf%catalog,1), ubound(fwdModelConf%catalog,1), 2
      do c = 1, size(fwdModelConf%catalog,2)
        ! Note that we don't deallocate the signals/sidebands stuff for each line
        ! as these are shallow copies of the main spectroscopy catalog stuff
        call deallocate_test ( fwdModelConf%catalog(s,c)%lines, &
          & 'fwdModelConf%catalog(?,?)%lines', moduleName )
        call deallocate_test ( fwdModelConf%catalog(s,c)%polarized, &
          & 'fwdModelConf%catalog(?,?)%polarized', moduleName )
      end do
    end do

    deallocate ( fwdModelConf%catalog, stat=s )
    if ( s /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'fwdModelConf%catalog' )

  end subroutine Destroy_Species_Data

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Species_Data_m

! $Log$
! Revision 2.19  2004/11/04 03:40:54  vsnyder
! Index spectroscopy catalog by molecule instead of searching
!
! Revision 2.18  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
