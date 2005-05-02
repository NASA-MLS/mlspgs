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
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DumpCommand ( Root, QuantityTemplatesDB, &
    & VectorTemplates, Vectors, ForwardModelConfigs, HGrids, VGrids )

  ! Process a "dump" command

    use AntennaPatterns_m, only: Dump_Antenna_Patterns_Database
    use Calendar, only: Duration_Formatted, Time_T, TK
    use Declaration_table, only: Num_Value
    use Expr_m, only: Expr
    use FilterShapes_m, only: Dump_Filter_Shapes_Database, &
      & Dump_DACS_Filter_Database
    use ForwardModelConfig, only: Dump, ForwardModelConfig_T
    use HGridsDatabase, only: Dump, HGRID_T
    use Init_Tables_Module, only: F_AllForwardModels, F_AllHGrids, F_AllLines, &
      & F_AllPFA, F_AllQuantityTemplates, F_AllSignals, F_AllSpectra, &
      & F_AllVectors, F_AllVectorTemplates, F_AllVGrids, F_AntennaPatterns, &
      & F_Details, F_DACSFilterShapes, F_FilterShapes, F_ForwardModel, F_HGrid, &
      & F_Lines, F_Mark, F_PfaData, F_PfaFiles, F_PointingGrids, F_Quantity, &
      & F_Signals,  F_Spectroscopy, F_Stop, F_Template, F_Text, F_TGrid, &
      & F_Vector, F_VGrid, S_Quantity, S_VectorTemplate
    use Intrinsic, only: PHYQ_Dimensionless
    use MoreTree, only: Get_Boolean, Get_Field_ID, Get_Spec_ID
    use MLSSignals_m, only: Dump, Signals
    use Output_M, only: Output
    use PFADataBase_m, only: Dump, Dump_PFADataBase, Dump_PFAFileDataBase, PFAData
    use PointingGrid_m, only: Dump_Pointing_Grid_Database
    use QuantityTemplates, only: Dump, QuantityTemplate_T
    use SpectroscopyCatalog_m, only: Catalog, Dump, Dump_Lines_Database, Lines
    use String_Table, only: Get_String
    use Tree, only: Decoration, Node_Id, Nsons, Source_Ref, Sub_rosa, Subtree
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

    real(tk) :: CPUTime, CPUTimeBase = 0.0_tk
    character(8) :: Date
    integer :: Details
    integer :: FieldIndex
    logical :: GotOne ! of something -- used to test loop completion
    integer :: GSON, I, J, K, L, Look
    integer :: QuantityIndex
    integer :: Son
    integer :: Source ! column*256 + line
    character :: TempText*20, Text*255
    type(time_t) :: Time
    character(10) :: TimeOfDay
    integer :: VectorIndex
    integer :: Type     ! of the Details expr -- has to be num_value
    integer :: Units(2) ! of the Details expr -- has to be phyq_dimensionless
    double precision :: Values(2) ! of the Details expr
    integer :: What

    ! Error codes
    integer, parameter :: Dimless = 1
    integer, parameter :: NoFWM = dimless + 1
    integer, parameter :: NoHGrid = NoFWM + 1
    integer, parameter :: NoLines = noHGrid + 1
    integer, parameter :: NoQT = noLines + 1
    integer, parameter :: NoSignals = noQT + 1
    integer, parameter :: NoTG = noSignals + 1
    integer, parameter :: NoVectors = noTG + 1
    integer, parameter :: NoVG = noVectors + 1
    integer, parameter :: NoVT = noVG + 1
    integer, parameter :: Numeric = noVT + 1
    integer, parameter :: Stop = numeric + 1
    integer, parameter :: Unknown = stop + 1 ! Unknown template

    details = 0
    do j = 2, nsons(root)
      son = subtree(j,root) ! The argument
      fieldIndex = get_field_id(son)
      if (nsons(son) > 1) gson = subtree(2,son) ! Now value of said argument
      select case ( fieldIndex )
      case ( f_allForwardModels, f_allHGrids, f_allLines, f_allPFA, &
        & f_allQuantityTemplates, f_allSignals, f_allSpectra, f_allVectors, &
        & f_allVectorTemplates, f_allVGrids, f_antennaPatterns, &
        & f_DACSfilterShapes, f_filterShapes, f_pfaFiles, f_pointingGrids, &
        & f_stop )
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
          case ( f_allLines )
            if ( associated(lines) ) then
              call dump_lines_database
            else
              call announceError ( son, noLines )
            end if
          case ( f_allPFA )
            call Dump_PFADataBase ( details )
          case ( f_allQuantityTemplates )
            if ( present(quantityTemplatesDB) ) then
              call dump ( quantityTemplatesDB )
            else
              call announceError ( son, noQT )
            end if
          case ( f_allSignals )
            if ( associated(signals) ) then
              call dump ( signals, details=details>0 )
            else
              call announceError ( son, noSignals )
            end if
          case ( f_allSpectra )
            call dump ( catalog, details=details )
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
          case ( f_pfaFiles )
            call dump_PFAFileDatabase ( details )
          case ( f_pointingGrids )
            call dump_pointing_grid_database ( son )
          case ( f_stop )
            call announceError ( son, stop )
            stop
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
      case ( f_lines )
        do i = 2, nsons(son)
          what = decoration(decoration(subtree(i,son)))
          call output ( what, after=': ' )
          call dump ( lines(what) )
        end do
      case ( f_mark )
        if ( get_boolean(son) ) call cpu_time ( cpuTimeBase )
      case ( f_pfaData )
        do i = 2, nsons(son)
          look = decoration(decoration(subtree(i,son)))
          call dump ( pfaData(look), details, look )
        end do
      case ( f_quantity ) ! Dump vector quantities
        do i = 2, nsons(son)
          gson = subtree(i,son)
          vectorIndex = decoration(decoration(subtree(1,gson)))
          quantityIndex = decoration(decoration(decoration(subtree(2,gson))))
          call dump ( GetVectorQtyByTemplateIndex( &
            & vectors(vectorIndex), quantityIndex), details=details )
        end do
      case ( f_signals )
        do i = 2, nsons(son)
          what = decoration(decoration(subtree(i,son)))
          call output ( what, after=': ' )
          call dump ( signals(what), details=details>0 )
        end do
      case ( f_spectroscopy )
        do i = 2, nsons(son)
          what = decoration(subtree(i,son))
          call output ( what, after=': ' )
          call dump ( catalog(what), details=details )
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
      case ( f_text )
        do k = 2, nsons(son)
          call get_string ( sub_rosa(subtree(k,son)), text, strip=.true. )
          ! Replace format marks: %[Cc] => CPU time in YyDdHH:MM:SS.SSS format
          !                       %[Dd] => date and time of day
          !                       %[Ss] => CPU time in seconds
          !                       %[Tt] => time of day
          gotOne = .false.
          do
            i = max(index(text,'%c'), index(text,'%C'))
            if ( i /= 0 ) then
              tempText = ''
              l = 1 ! Position in TempText
              call cpu_time ( cpuTime )
              time = duration_formatted ( cpuTime - cpuTimeBase )
              if ( time%year /= 0 ) then
                gotOne = .true.
                write ( tempText, * ) time%year
                tempText = adjustl(tempText)
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = 'y'
              end if
              if ( time%day /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), * ) time%day
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = 'd'
              end if
              if ( time%hours /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), '(i2.2)' ) time%hours
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = ':'
              end if
              if ( time%minutes /= 0 .or. gotOne ) then
                gotOne = .true.
                write ( tempText(l:), '(i2.2)' ) time%minutes
                tempText(l:) = adjustl(tempText(l:))
                l = len_trim(tempText) + 2
                tempText(l-1:l-1) = ':'
              end if
              write ( tempText(l:), '(f6.3)' ) time%seconds
              tempText(l:) = adjustl(tempText(l:))
              l = len_trim(tempText)
              text = text(:i-1) // tempText(:l) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%d'), index(text,'%D'))
            if ( i /= 0 ) then
              call date_and_time ( date=date, time=timeOfDay )
              text = text(:i-1) // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // &
                & ' ' // timeOfDay(1:2) // ':' // timeOfDay(3:4) // ':' // timeOfDay(5:10) // &
                & text(i+2:)
              cycle
            end if
            i = max(index(text,'%l'), index(text,'%L'))
            if ( i /= 0 ) then
              source = source_ref(subtree(k,son))
              write ( tempText(1:10), * ) source/256
              write ( tempText(11:20), * ) mod(source,256)
              text = text(:i-1) // "line " // trim(adjustl(tempText(1:10))) // &
                & ", column " // trim(adjustl(tempText(11:20))) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%s'), index(text,'%S'))
            if ( i /= 0 ) then
              write ( tempText, * ) cpuTime - cpuTimeBase
              text = text(:i-1) // trim(adjustl(tempText)) // text(i+2:)
              cycle
            end if
            i = max(index(text,'%t'), index(text,'%T'))
            if ( i /= 0 ) then
              call date_and_time ( time=timeOfDay )
              text = text(:i-1) // timeOfDay(1:2) // ':' // timeOfDay(3:4) // &
                & ':' // timeOfDay(5:10) // text(i+2:)
              cycle
            end if
            exit ! Didn't find a format trigger
          end do
          call output ( trim(text), advance='yes' )
        end do ! k
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
      use Output_m, only: NewLine

      integer, intent(in) :: What, Where

      call StartErrorMessage ( where )

      select case ( what )
      case ( dimless )
        call output ( "The details field is not unitless." )
      case ( noFWM )
        call output ( "Can't dump Forward Model Configs here." )
      case ( noHGrid )
        call output ( "Can't dump HGrids here." )
      case ( noLines )
        call output ( "Can't dump Lines here." )
      case ( noQT )
        call output ( "Can't dump Quantity Templates here." )
      case ( noSignals )
        call output ( "Can't dump Signals here." )
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
      case ( stop )
        call output ( "Program stopped by /stop field on DUMP statement." )
      case ( unknown )
        call output ( "Can't figure out what kind of template it is." )
      end select
      call newLine
    end subroutine AnnounceError

  end subroutine DumpCommand

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DumpCommand_M

! $Log$
! Revision 2.22  2005/05/02 23:11:37  vsnyder
! Add dump of PFAFiles database
!
! Revision 2.21  2005/04/01 20:48:28  vsnyder
! Add mark and text fields to dump command
!
! Revision 2.20  2005/03/26 01:34:00  vsnyder
! Add stop message
!
! Revision 2.19  2005/03/15 01:36:08  vsnyder
! Add newline after error messages
!
! Revision 2.18  2005/01/12 03:18:51  vsnyder
! Add item number to PFA dump
!
! Revision 2.17  2004/12/28 00:22:03  vsnyder
! Add not_used_here
!
! Revision 2.16  2004/12/13 20:13:04  vsnyder
! Add dumps for AllLines, AllSignals, AllSpectra, Lines, Signals, Spectroscopy,
! and a Stop command.
!
! Revision 2.15  2004/11/04 06:37:34  vsnyder
! Index spetroscopy catalog by molecule instead of searching
!
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
