program mockup
   use CFM_MLSSetup_m, only: CFM_MLSSetup, CFM_MLSCleanup
   use CFM_VGrid, only: CreateVGrid
   use Chunks_m, only: MLSChunk_T
   use ForwardModelConfig, only: ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T
   use VGridsDatabase, only: VGrid_T, Dump
   use Intrinsic, only: phyq_pressure, l_zeta
   use Init_Tables_Module, only: l_logarithmic

   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

   integer :: error
   type(ForwardModelConfig_T), dimension(:), pointer :: forwardModelConfigDatabase
   type(MLSFile_T), dimension(:), pointer :: filedatabase
   type(MLSChunk_T) :: fakeChunk
   type(VGrid_T) :: vGridStandard

   call CFM_MLSSetup(error, ForwardModelConfigDatabase, filedatabase, fakeChunk)
   if (error /=0) stop

   vGridStandard = CreateVGrid(l_zeta, l_logarithmic, &
                               1000.0d0, "37:6", phyq_pressure)

   call dump(vGridStandard)

   call CFM_MLSCleanup
end program
