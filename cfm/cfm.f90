! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module cfm          ! callable forward model
  ! This module merely unifies all the .mod files that
  ! will be USEd by any fortran90 program that is compiled
  ! and with the callable forward model and links against libcfm.a or libcfm_all.a
   use CFM_MLSSetup_m, only: CFM_MLSCleanup, CFM_MLSSetup
   use Chunks_m, only: MLSChunk_T
   use CFM_IO_M, only: Read_Spectroscopy, ReadDACSFilterShapes, &
                     ReadAntennaPatterns, ReadFilterShapes, &
                     ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC
   use CFM_Matrix_M, only: CreatePlainMatrix
   use PFADatabase_m, only: Destroy_PFADataBase
   use L2PC_m, only: DestroyL2PCDatabase
   use FilterShapes_m, only: Destroy_DACS_Filter_Database, &
                             Destroy_Filter_Shapes_Database
   use AntennaPatterns_m, only: Destroy_Ant_Patterns_Database
   use SpectroscopyCatalog_m, only: Destroy_SpectCat_Database, &
                                    Destroy_Line_Database
   use PointingGrid_m, only: Destroy_Pointing_Grid_Database
   use CFM_VGrid_m, only: CreateVGrid
   use VGridsDatabase, only: DestroyVGridContents, VGrid_T, Dump
   use CFM_HGrid_m, only: CreateRegularHGrid
   use HGridsDatabase, only: DestroyHGridContents, HGrid_T, Dump
   use CFM_FGrid_m, only: CreateFGrid
   use FGrid, only: FGrid_T, DestroyFGridContents, Dump
   use CFM_QuantityTemplate_m, only: CreateQtyTemplate
   use QuantityTemplates, only: Dump, AddQuantityTemplateToDatabase, &
                                DestroyQuantityTemplateDatabase, &
                                QuantityTemplate_T
   use CFM_VectorTemplate_m, only: CreateVectorTemplate
   use VectorsModule, only: Dump, Vector_T, VectorValue_T, &
                            VectorTemplate_T, DestroyVectorTemplateInfo, &
                            DestroyVectorInfo, GetVectorQtyByTemplateIndex, &
                            operator(+), operator(-)
   use CFM_Vector_m, only: CreateVector
   use CFM_Fill_m, only: ExplicitFillVectorQuantity, ApplyBaseline, &
                         FillVectorQuantityFromL1B, FillPhitanQuantity, &
                         SpreadFillVectorQuantity, FillPtanQuantity
   use ForwardModelConfig, only: ForwardModelConfig_T, Dump
   use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
   use CFM_FWDMDL_M, only: ForwardModel
   use MLSCommon, only: MLSFile_T, r8
   use Init_tables_module, only: l_logarithmic, l_zeta, l_temperature, &
                                 L_IntermediateFrequency, l_vmr, l_gph, &
                                 l_ptan, l_radiance, l_orbitInclination, &
                                 l_tngtgeodalt, l_tngtgeocalt, l_o3, &
                                 phyq_pressure, phyq_angle, l_h2o, l_refgph, &
                                 l_phitan, l_explicit, l_l1bMAFBaseline, &
                                 l_frequency, l_usbFrequency, l_lsbFrequency, &
                                 l_ghz
   use MLSFiles, only: GetMLSFileByType, InitializeMLSFile, mls_openFile, &
                       AddFileToDatabase
   use Intrinsic, only: l_hdf
   use Hdf, only: DFACC_RDONLY
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error
   use MatrixModule_1, only: Matrix_T

   implicit none
   public

contains
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
    character (len=*), parameter :: ModuleName= &
         "$RCSfile$"
    character (len=*), parameter :: IdParm = &
         "$Id$"
    character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module cfm

! $Log$
! Revision 1.3  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.2  2010/06/29 15:29:34  honghanh
! Develop FillPtanQuantity to compute ptan, instead of using
! Get2DHydrostaticTangentPressure
!
! Revision 1.1  2010/06/29 02:34:03  honghanh
! Change cfm_m.f90 to cfm.f90 to match the naming scheme of cfm directory
!
! Revision 1.4  2010/06/29 02:28:17  honghanh
! Change mockup to import functions and literals from CFM module
!
! Revision 1.2  2010/06/16 20:23:52  honghanh
! Update cfm_fwdmdl to cfm_fwdmdl_m
!
! Revision 1.1  2010/06/03 23:31:57  pwagner
! First commit
!
