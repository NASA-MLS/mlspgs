! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DumpCommand_M

! Process a "dump" command.

  implicit NONE
  private

  public :: DumpCommand


!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine DumpCommand ( Root, QuantityTemplatesDB, &
    & VectorTemplates, Vectors, ForwardModelConfigs, HGrids, VGrids )

  ! Process a "dump" command

  ! The fields can be Quantity to dump a vector quantity, Vector to dump an
  ! entire vector, Template to dump a quantity template or vector template,
  ! ForwardModel to dump a forward model config, or Details to specify
  ! the level of detail for subsequent dumps.

    use AntennaPatterns_m, only: Dump_Antenna_Patterns_Database
    use Declaration_table, only: Num_Value
    use Expr_m, only: Expr
    use FilterShapes_m, only: Dump_Filter_Shapes_Database, &
      & Dump_DACS_Filter_Database
    use ForwardModelConfig, only: Dump, ForwardModelConfig_T
    use HGridsDatabase, only: Dump, HGRID_T
    use Init_Tables_Module, only: F_AllForwardModels, F_AllHGrids, F_AllPFA, &
      & F_AllQuantityTemplates, F_AllVectors, F_AllVectorTemplates, &
      & F_AllVGrids, F_AntennaPatterns, F_Details, F_DACSFilterShapes, &
      & F_FilterShapes, F_ForwardModel, F_HGrid, F_PfaData, &
      & F_PointingGrids, F_Quantity, F_Spectroscopy, F_Template, F_TGrid, &
      & F_Vector, F_VGrid, S_Quantity, S_VectorTemplate
    use Intrinsic, only: PHYQ_Dimensionless
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
    use Output_M, only: Output
    use PFADataBase_m, only: Dump, Dump_PFADataBase, PFAData
    use PointingGrid_m, only: Dump_Pointing_Grid_Database
    use QuantityTemplates, only: Dump, QuantityTemplate_T
    use SpectroscopyCatalog_m, only: Catalog, Dump
    use Tree, only: Decoration, Node_Id, Nsons, Subtree
    use Tree_Types, only: N_Spec_Args
    use VectorsModule, only: Dump, & ! for vectors, vector quantities and templates
      & GetVectorQtyByTemplateIndex, Vector_T, VectorTemplate_T
    use VGridsDatabase, only: Dump, VGrid_T

    integer, intent(in) :: Root ! Root of the parse tree for the dump command
    ! Databases:
    type (quantityTemplate_t), dimension(:), intent(in), optional :: QuantityTemplatesDB
    type (vectorTemplate_T), dimension(:), intent(in), optional :: VectorTemplates
    type (vector_T), dimension(:), intent(in), optional :: Vectors
    type (forwardModelConfig_t), dimension(:), pointer, optional :: ForwardModelConfigs
    type (HGrid_T), dimension(:), intent(in), optional :: HGrids
    type (VGrid_T), dimension(:), pointer, optional :: VGrids

    integer :: Details
    integer :: FieldIndex
    integer :: GSON, I, J, Look
    integer :: QuantityIndex
    integer :: Son
    integer :: VectorIndex
    integer :: Type     ! of the Details expr -- has to be num_value
    integer :: Units(2) ! of the Details expr -- has to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr
    integer :: What

    ! Error codes
    integer, parameter :: Dimless = 1
    integer, parameter :: NoFWM = dimless + 1
    integer, parameter :: NoHGrid = NoFWM + 1
    integer, parameter :: NoQT = noHGrid + 1
    integer, parameter :: NoSpec = noQT + 1
    integer, parameter :: NoTG = NoSpec + 1
    integer, parameter :: NoVectors = noTG + 1
    integer, parameter :: NoVG = noVectors + 1
    integer, parameter :: NoVT = noVG + 1
    integer, parameter :: Numeric = noVT + 1
    integer, parameter :: Unknown = numeric + 1 ! Unknown template

    details = 0
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case ( f_allForwardModels, f_allHGrids, f_allPFA, f_allQuantityTemplates, &
        & f_allVectors, f_allVectorTemplates, f_allVGrids, f_antennaPatterns, &
        & f_DACSfilterShapes, f_filterShapes, f_pointingGrids, f_spectroscopy )
        if ( get_boolean(son) ) then
          select case ( fieldIndex )
          case ( f_allForwardModels )
            if ( present(forwardModelConfigs) ) then
              call dump ( forwardModelConfigs, where=son )
            else
              call announceError ( son, noFWM )
            end if
          case ( f_allHGrids )
            if ( present(hGrids) ) then
              call dump ( hGrids )
            else
              call announceError ( son, noHGrid )
            end if
          case ( f_allPFA )
            call Dump_PFADataBase ( details )
          case ( f_allQuantityTemplates )
            if ( present(quantityTemplatesDB) ) then
              call dump ( quantityTemplatesDB )
            else
              call announceError ( son, noQT )
            end if
          case ( f_allVectors )
            if ( present(vectors) ) then
              call dump ( vectors )
            else
              call announceError ( son, noVectors )
            end if
          case ( f_allVectorTemplates )
            if ( present(vectorTemplates) ) then
              call dump ( vectorTemplates )
            else
              call announceError ( son, noVT )
            end if
          case ( f_allVGrids )
            if ( present(vGrids) ) then
              call dump ( vGrids )
            else
              call announceError ( son, noVG )
            end if
          case ( f_antennaPatterns )
            call dump_antenna_patterns_database ( son )
          case ( f_DACSfilterShapes )
            call dump_dacs_filter_database ( son )
          case ( f_filterShapes )
            call dump_filter_shapes_database ( son )
          case ( f_pointingGrids )
            call dump_pointing_grid_database ( son )
          case ( f_spectroscopy )
            if ( associated(catalog) ) then
              call dump ( catalog, details=details )
            else
              call announceError ( son, noSpec )
            end if
          end select
        end if
      case ( f_details )
        call expr ( gson, units, values, type )
        if ( units(1) /= phyq_dimensionless ) call AnnounceError ( gson, dimless )
        if ( type /= num_value ) call announceError ( gson, numeric )
        details = nint(values(1))
      case ( f_forwardModel ) ! Dump forward model configs
        if ( present(forwardModelConfigs) ) then
          do i = 2, nsons(son)
            call dump ( & ! has no details switch
              & forwardModelConfigs(decoration(decoration(subtree(i,son)))) )
          end do
        else
          call announceError ( gson, noFWM )
        end if
      case ( f_hGrid )    ! Dump HGrids
        if ( present(hGrids) ) then
          do i = 2, nsons(son)
            call output ( ' HGrid ' )
            call dump ( & ! has no details switch
              & hGrids(decoration(decoration(subtree(i,son)))) )
          end do
        else
          call announceError ( gson, noHGrid )
        end if
      case ( f_pfaData )
        do i = 2, nsons(son)
          call dump ( pfaData(decoration(decoration(subtree(i,son)))), details )
        end do
      case ( f_quantity ) ! Dump vector quantities
        do i = 2, nsons(son)
          gson = subtree(i,son)
          vectorIndex = decoration(decoration(subtree(1,gson)))
          quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          call dump ( GetVectorQtyByTemplateIndex( &
            & vectors(vectorIndex), quantityIndex), details=details )
        end do
      case ( f_template ) ! Dump vector templates or quantity templates
        do i = 2, nsons(son)
          gson = subtree(i,son)
          look = decoration(gson)
          if ( node_id(look) /= n_spec_args ) call announceError ( gson, unknown )
          what = decoration(look)
          select case ( get_spec_id(look) )
          case ( s_quantity )
            if ( present(quantityTemplatesDB) ) then
              call output ( ' Quantity template' )
              call dump ( quantityTemplatesDB(what), details=details )
            else
              call announceError ( gson, noQT )
            end if
          case ( s_vectorTemplate )
            if ( present(vectorTemplates) ) then
              call output ( ' Vector template' )
              call dump ( vectorTemplates(what), details=details, quantities=quantityTemplatesDB )
            else
              call announceError ( gson, noVT )
            end if
          end select
        end do
      case ( f_tGrid )
        if ( present(vGrids) ) then
          do i = 2, nsons(son)
            call output ( ' TGrid ' )
            call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
          end do
        else
          call announceError ( gson, noTG )
        end if
      case ( f_vector ) ! Dump entire vectors
        if ( present(vectors) ) then
          do i = 2, nsons(son)
            call output ( ' Vector ' )
            call dump ( vectors(decoration(decoration(subtree(i,son)))), details=details )
          end do
        else
          call announceError ( gson, noVectors )
        end if
      case ( f_vGrid )
        if ( present(vGrids) ) then
          do i = 2, nsons(son)
            call output ( ' VGrid ' )
            call dump ( vGRids(decoration(decoration(subtree(i,son)))), details=details )
          end do
        else
          call announceError ( gson, noVG )
        end if
      end select
    end do

  contains

    subroutine AnnounceError ( where, what )
      use MoreTree, only: StartErrorMessage

      integer, intent(in) :: What, Where

      call StartErrorMessage ( where )

      select case ( what )
      case ( dimless )
        call output ( "The details field is not unitless." )
      case ( noFWM )
        call output ( "Can't dump Forward Model Configs here." )
      case ( noHGrid )
        call output ( "Can't dump HGrids here." )
      case ( noQT )
        call output ( "Can't dump Quantity Templates here." )
      case ( noSpec )
        call output ( "Can't dump spectroscopy database here." )
      case ( noTG )
        call output ( "Can't dump TGrids here." )
      case ( noVectors )
        call output ( "Can't dump Vectors here." )
      case ( noVG )
        call output ( "Can't dump VGrids here." )
      case ( noVT )
        call output ( "Can't dump Vector Templates here." )
      case ( numeric )
        call output ( "The details field is not numeric." )
      case ( unknown )
        call output ( "Can't figure out what kind of template it is." )
      end select
    end subroutine AnnounceError

  end subroutine DumpCommand

end module DumpCommand_M

! $Log$
! Revision 2.14  2004/11/01 20:16:20  vsnyder
! Check for spectroscopy catalog before trying to dump it
!
! Revision 2.13  2004/10/30 00:26:46  vsnyder
! Add 'spectroscopy' field to DumpCommand
!
! Revision 2.12  2004/10/06 20:19:39  vsnyder
! Cannonball polishing
!
! Revision 2.11  2004/09/24 22:24:20  vsnyder
! Make PFA dump aware of 'details' switch
!
! Revision 2.10  2004/07/17 02:28:19  vsnyder
! Add dump for entire PFA database
!
! Revision 2.9  2004/06/12 00:41:30  vsnyder
! Allow all fields except details to be arrays
!
! Revision 2.8  2004/06/09 19:59:38  pwagner
! Gets PFAData type and dump method from PFADataBase_m
!
! Revision 2.7  2004/06/08 20:20:18  vsnyder
! Add tGrid
!
! Revision 2.6  2004/05/29 02:50:49  vsnyder
! Added more dumps
!
! Revision 2.5  2004/05/22 02:31:23  vsnyder
! Dump PFAData, VGrids
!
! Revision 2.4  2004/05/20 19:47:36  vsnyder
! Move Dump*Hgrid from Dumper to HgridsDatabse
!
! Revision 2.3  2004/05/18 01:18:51  vsnyder
! Add dump for HGrid
!
! Revision 2.2  2004/05/11 02:53:29  vsnyder
! Remove USEs for unreferenced symbols
!
! Revision 2.1  2004/05/01 04:04:16  vsnyder
! Initial commit
!
