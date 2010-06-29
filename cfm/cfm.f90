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
   use CFM_MLSSetup_m, only: CFM_MLSCleanup, MLSChunk_T, CFM_MLSSetup
   use CFM_IO_M, only: Read_Spectroscopy, ReadDACSFilterShapes, &
                     ReadAntennaPatterns, ReadFilterShapes, &
                     ReadPointingGrids, ReadPFAFile, ReadHDF5L2PC, &
                     Destroy_DACS_Filter_Database, &
                     Destroy_Filter_Shapes_Database, &
                     Destroy_Ant_Patterns_Database, &
                     Destroy_SpectCat_Database, &
                     Destroy_Line_Database, &
                     Destroy_Pointing_Grid_Database, &
                     DestroyL2PCDatabase, Destroy_PFADataBase
   use CFM_VGrid_m, only: CreateVGrid, DestroyVGridContents, &
                        VGrid_T, Dump
   use CFM_HGrid_m, only: CreateRegularHGrid, HGrid_T, &
                        DestroyHGridContents, Dump
   use CFM_FGrid_m, only: CreateFGrid, FGrid_T, DestroyFGridContents, &
                        Dump
   use CFM_QuantityTemplate_m, only: CreateQtyTemplate, Dump, &
                        AddQuantityTemplateToDatabase, &
                        DestroyQuantityTemplateDatabase, &
                        QuantityTemplate_T
   use CFM_VectorTemplate_m, only: CreateVectorTemplate, Dump, &
                        VectorTemplate_T, DestroyVectorTemplateInfo
   use CFM_Vector_m, only: CreateVector, Dump, &
                         Vector_T, VectorValue_T, &
                         DestroyVectorInfo, GetVectorQtyByTemplateIndex
   use CFM_Fill_m, only: ExplicitFillVectorQuantity, &
                         FillVectorQuantityFromL1B, FillPhitanQuantity, &
                         SpreadFillVectorQuantity
   use CFM_FWDMDL_M, only: ForwardModel, FORWARDMODELSTATUS_T, &
                         ForwardModelConfig_T
   use MLSCommon, only: MLSFile_T, r8
   use Init_tables_module, only: l_logarithmic, l_zeta, l_temperature, &
                                 L_IntermediateFrequency, l_vmr, l_gph, &
                                 l_ptan, l_radiance, l_orbitInclination, &
                                 l_tngtgeodalt, l_tngtgeocalt, l_o3, &
                                 phyq_pressure, phyq_angle, l_h2o, l_refgph, &
                                 l_phitan, l_explicit
   use MLSFiles, only: GetMLSFileByType, InitializeMLSFile, mls_openFile, &
                       AddFileToDatabase
   use ScanModelModule, only: Get2DHydrostaticTangentPressure
   use Intrinsic, only: l_hdf
   use Hdf, only: DFACC_RDONLY
   use MLSMessageModule, only: MLSMessage, MLSMSG_Error

   implicit none
   public

!---------------------------- RCS Ident Info ------------------------------
   character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
   character (len=*), parameter :: IdParm = &
       "$Id$"
   character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
contains
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module cfm

! $Log$
! Revision 1.4  2010/06/29 02:28:17  honghanh
! Change mockup to import functions and literals from CFM module
!
! Revision 1.2  2010/06/16 20:23:52  honghanh
! Update cfm_fwdmdl to cfm_fwdmdl_m
!
! Revision 1.1  2010/06/03 23:31:57  pwagner
! First commit
!
