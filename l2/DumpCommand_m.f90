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
  character (len=len(idParm)), private :: Id = idParm
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

    use Declaration_table, only: Num_Value
    use Expr_m, only: Expr
    use ForwardModelConfig, only: Dump, ForwardModelConfig_T
    use HGridsDatabase, only: Dump, HGRID_T
    use Init_Tables_Module, only: F_Details, F_ForwardModel, F_HGrid, &
      & F_PfaData, F_Quantity, F_Template, F_Vector, F_VGrid, &
      & S_Quantity, S_VectorTemplate
    use Intrinsic, only: PHYQ_Dimensionless
    use MoreTree, only: Get_Field_ID, Get_Spec_ID
    use Output_M, only: Output
    use PFAData_m, only: Dump, PFAData
    use QuantityTemplates, only: Dump, QuantityTemplate_T
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
    type (forwardModelConfig_t), dimension(:), intent(in), optional :: ForwardModelConfigs
    type (HGrid_T), dimension(:), intent(in), optional :: HGrids
    type (VGrid_T), dimension(:), intent(in), optional :: VGrids

    integer :: Details
    integer :: FieldIndex
    integer :: GSON, J, Look
    integer :: QuantityIndex
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
    integer, parameter :: NoVectors = noQT + 1
    integer, parameter :: NoVG = noVectors + 1
    integer, parameter :: NoVT = noVG + 1
    integer, parameter :: Numeric = noVT + 1
    integer, parameter :: Unknown = numeric + 1 ! Unknown template

    details = 0
    do j = 2, nsons(root)
      gson = subtree(j,root) ! The argument
      fieldIndex = get_field_id(gson)
      if (nsons(gson) > 1) gson = subtree(2,gson) ! Now value of said argument
      select case ( fieldIndex )
      case ( f_details )
        call expr ( gson, units, values, type )
        if ( units(1) /= phyq_dimensionless ) call AnnounceError ( gson, dimless )
        if ( type /= num_value ) call AnnounceError ( gson, numeric )
        details = nint(values(1))
      case ( f_forwardModel ) ! Dump a forward model config
        if ( present(forwardModelConfigs) ) then
          what = decoration(decoration(gson))
          call dump ( forwardModelConfigs(what) ) ! has no details switch
        else
          call announceError ( gson, noFWM )
        end if
      case ( f_hGrid )    ! Dump an HGrid
        if ( present(hGrids) ) then
          what = decoration(decoration(gson))
          call output ( ' HGrid ' )
          call dump ( hGrids(what) ) ! has no details switch
        else
          call announceError ( gson, noHGrid )
        end if
      case ( f_pfaData )
        call output ( ' PFA Data', advance='yes' )
        call dump ( pfaData(decoration(decoration(gson))) )
      case ( f_quantity ) ! Dump a vector quantity
        vectorIndex = decoration(decoration(subtree(1,gson)))
        quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
        call dump ( GetVectorQtyByTemplateIndex( &
          & vectors(vectorIndex), quantityIndex), details=details )
      case ( f_template ) ! Dump a vector template or a quantity template
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
      case ( f_vector ) ! Dump an entire vector
        if ( present(vectors) ) then
          call output ( ' Vector ' )
          call dump ( vectors(decoration(decoration(gson))), details=details )
        else
          call announceError ( gson, noVectors )
        end if
      case ( f_vGrid )
        if ( present(vGrids) ) then
          call output ( ' VGrid ' )
          call dump ( vGRids(decoration(decoration(gson))), details=details )
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
