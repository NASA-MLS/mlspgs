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
   use CFM_MLSSetup_m, only: CFM_MLSCleanup, CFM_MLSSetup, CreateMLSValue_O2, &
                             CreateMLSValue_EarthReflectivity, &
                             CreateMLSValue_LSF, CreateMLSValue_FromL1BOA, &
                             CreateMLSValue_SpaceRadiance, GetConstantQuantities, &
                             CreateMLSValue_ElevationOffset, &
                             timeRange2MafRange
   use Chunks_m, only: MLSChunk_T
   use CFM_IO_M, only: Read_Spectroscopy, ReadDACSFilterShapes, &
                     ReadAntennaPatterns, ReadFilterShapes, &
                     ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC, &
                     Read_ptan, Write_To_File1, Write_To_File2, & ! Added by Pranjit Saha
                     read_H5inputdata, & ! Added by Zheng Qu
                     read_txtinputdata, & ! Added by Zheng Qu
                     Write_To_HDF5, Write2HDF5, &  ! Added by Zheng Qu
                     Copy_RadianceJacobian, path_join,  & ! Added by Zheng Qu
                     Filter_RadianceJacobian, log_inputdata, & ! Added by Zheng Qu
                     runLogFileUnit, errorLogFileUnit, runlog, & ! Added by Zheng Qu
                        errlog, openLogFiles  ! Added by Zheng Qu
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
                                QuantityTemplate_T, DestroyQuantityTemplateContents
   use CFM_VectorTemplate_m, only: CreateVectorTemplate
   use VectorsModule, only: Dump, Vector_T, VectorValue_T, &
                            VectorTemplate_T, DestroyVectorTemplateInfo, &
                            DestroyVectorInfo, GetVectorQtyByTemplateIndex, &
                            operator(+), operator(-)
   use CFM_Vector_m, only: CreateVector, CreateValue4AgileVector, &
                           CreateAgileVector, AddValue2Vector, &
                           DestroyAgileVectorContent, DestroyVectorValueContent, &
                           CloneAgileVector
   use CFM_Fill_m, only: ExplicitFillVectorQuantity, ApplyBaseline, &
                         FillVectorQuantityFromL1B, FillPhitanQuantity, &
                         SpreadFillVectorQuantity, FillPtanQuantity
   use ForwardModelConfig, only: ForwardModelConfig_T, Dump
   use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
   use CFM_FWDMDL_M, only: ForwardModel, ForwardModel2
   use MLSCommon, only: MLSFile_T, r8
   use Init_tables_module ! essentially, all the l_..., and some phyq
   use MLSFiles, only: GetMLSFileByType, InitializeMLSFile, mls_openFile, &
                       AddFileToDatabase, mls_closeFile
   use Intrinsic, only: l_hdf
   use Hdf, only: DFACC_RDONLY
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
                               MLSMSG_Debug, MLSMSG_Info, MLSMSG_Crash, &
                               MLSMSG_Success, MLSMessageSetup, & 
                               MLSMessageClose , MLSMessageExit

   use MatrixModule_1, only: Matrix_T, Dump, DestroyMatrix

   implicit none
   public
   private :: not_used_here

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
! Revision 1.25  2011/12/25 03:39:06  honghanh
! Fix cloning and memory cleanup method
! in forwardmodel2 subroutines.
!
! Revision 1.24  2011/12/15 16:53:24  honghanh
! Correct the name of CreateMLSValue_EarthReflectivity
!
! Revision 1.23  2011/12/14 22:54:18  honghanh
! Add timeRange2MafRange method in CFM.
!
! Revision 1.22  2011/11/01 22:16:11  honghanh
! Add API to destroy individual QuantityTemplate_T and VectorValue_T
!
! Revision 1.21  2011/10/31 20:04:59  honghanh
! Change CreateMLSValue_FromL1B to CreateMLSValue_FromL1BOA.
!
! Revision 1.20  2011/10/20 00:23:07  honghanh
! Add elevation offset creation subroutine to public API
!
! Revision 1.19  2011/10/19 19:33:10  honghanh
! Adding DestroyAgileVectorContent to CFM_Vector_m
!
! Revision 1.16  2011/08/25 20:10:32  honghanh
! Add other log levels other than Error and Warning from MLSMessageModule
!
! Revision 1.15  2011/06/21 19:04:48  honghanh
! Adding ForwardModel2
!
! Revision 1.14  2011/03/24 15:16:46  honghanh
! Add new interfaces for creating vector and vector values without going through quantity template databases
!
! Revision 1.13  2011/03/10 18:18:02  pwagner
! Everything from init_tables_module back
!
! Revision 1.12  2011/01/31 17:45:43  honghanh
! Make not_used_here a private function
!
! Revision 1.11  2010/12/10 00:48:28  pwagner       
! Everything from init_tables_module now available
!       
! Revision 1.10  2010/11/18 19:04:13  honghanh  
! New example to run forward model with a single maf
!
! Revision 1.9  2010/09/28 14:42:42  honghanh
! Add call to forwardModel with jacobian
!
! Revision 1.8  2010/09/21 15:06:49  honghanh
! Separate the running forward model example and reading observed radiance example
!
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
