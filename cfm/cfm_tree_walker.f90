! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module CFM_Tree_Walker_m

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions in the order indicated.

   implicit none
   private

   public :: Walk_Tree

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

   subroutine Walk_Tree ( Root, First_Section, ForwardModelConfigDatabase )

      use ForwardModelConfig, only: ForwardModelConfig_T
      use Global_Settings, only: Set_Global_Settings
      use Init_Tables_Module, only: Z_GLOBALSETTINGS, Z_MLSSIGNALS
      use MLSCommon, only: MLSFile_T
      use SpectroscopyCatalog_m, only: Spectroscopy
      use Time_M, only: Time_Now
      use Toggles, only: GEN, TOGGLE
      use Trace_m, only: TRACE_BEGIN, TRACE_END
      use Tree, only: DECORATION, NSONS, SUBTREE
      use MLSSignals_M, only: MLSSignals
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error

      integer, intent(in) ::       Root          ! Root of the abstract syntax tree
      integer, intent(in) ::       First_Section ! Index of son of root of first n_cf
      type (ForwardModelConfig_T), pointer, optional :: ForwardModelConfigDatabase(:)

      type (MLSFile_T), dimension(:), pointer :: File_Data_Base => NULL()
      integer :: I             ! Loop inductor, subtree index
      integer :: Section_Index ! Z_MLSSignals, Z_GlobalSettings
      integer :: Son           ! Tree index of son of a tree vertex
      real :: T1, T2           ! For timing

      logical, parameter :: TOOLKIT = .false.

      if ( toggle(gen) ) &
         & call trace_begin ( 'Walk_Tree', subtree(first_section,root) )
      call time_now ( t1 )
      nullify ( File_Data_Base )

      ! Loop over the tree
      do i = first_section, nsons(root)
         son = subtree(i,root)
         section_index = decoration(subtree(1,son))
         select case ( section_index )
         case ( z_GlobalSettings )
            if (present(ForwardModelConfigDatabase)) then
               call set_global_settings ( son, forwardModelConfigDatabase, File_Data_Base )
            else
               call MLSMessage(MLSMSG_Error, moduleName, &
               "Need ForwardModelConfigDatabase argument")
            end if
         case ( z_MLSSignals )
            call MLSSignals ( son )
         case default
            call print_unsupported_section_error (son, section_index)
         end select
      end do

      ! We only needed File_Data_Base to keep Spectroscopy and
      ! Set_Global_Settings happy
      if (associated(file_data_base)) deallocate ( File_Data_Base)

      call time_now ( t2 )
      !call output ( t2-t1, before='Timing to set up forward model = ', advance='yes' )

      if ( toggle(gen) ) call trace_end ( 'Walk_Tree' )

   end subroutine Walk_Tree

   subroutine print_unsupported_section_error (son, section_index)
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use String_Table, only: Display_String
      use Intrinsic, only: Section_Indices

      integer, intent(in) :: son, section_index

      call startErrorMessage ( son )
      call output ( 'Section ' )
      call display_string ( section_indices(section_index) )
      call output ( ' is not supported in callable forward model', advance='yes' )
      call MLSMessage ( MLSMSG_Error, moduleName, 'Unsupported CF section' )
   end subroutine

!--------------------------- end bloc --------------------------------------
   logical function not_used_here()
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
      not_used_here = (id(1:1) == ModuleName(1:1))
      print *, Id ! .mod files sometimes change if PRINT is added
   end function not_used_here
!---------------------------------------------------------------------------

end module

! $Log$
! Revision 1.3  2011/08/27 13:56:33  honghanh
! Removing Spectroscopy section from the allowed section
!
! Revision 1.2  2010/05/23 02:12:38  honghanh
! Add print_unsupported_section_error
!
! Revision 1.1  2010/02/17 16:41:33  honghanh
! An example of using CFM_MLSSetup and creating a VGrid
!
