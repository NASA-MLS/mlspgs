! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module cfm_eo_m  ! This module is to help cfm_mlssetup subroutine.
                 ! Quantities allocated in here should be deallocate by cfm_mlssetup_m
    use MLSCommon, only: r8

    implicit none

    private

    public :: CreateElevationOffsets, FillElevationOffsets

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

    real(r8), dimension(4), parameter :: velev32L = (/0.0200_r8, 0.0215_r8, 0.0186_r8, 0.0172_r8/)
    real(r8), dimension(4), parameter :: velev33L = (/-0.0002_r8, 0.0003_r8, 0.0003_r8, 0.0001_r8/)
    real(r8), dimension(4), parameter :: velev34L = (/0.0057_r8, 0.0057_r8, 0.0060_r8, 0.0063_r8/)
    real(r8), dimension(4), parameter :: velev33U = (/0.0001_r8, 0.0002_r8, -0.0000_r8, -0.0001_r8/)

    integer :: elev1L = 0
    integer :: elev2L = 0
    integer :: elev3L = 0
    integer :: elev4L = 0
    integer :: elev5L = 0
    integer :: elev6L = 0
    integer :: elev7L = 0
    integer :: elev8L = 0
    integer :: elev9L = 0
    integer :: elev10L = 0
    integer :: elev11L = 0
    integer :: elev12L = 0
    integer :: elev13L = 0
    integer :: elev14L = 0
    integer :: elev15L = 0
    integer :: elev16L = 0
    integer :: elev17L = 0
    integer :: elev18L = 0
    integer :: elev19L = 0
    integer :: elev20L = 0
    integer :: elev21L = 0
    integer :: elev22L = 0
    integer :: elev23L = 0
    integer :: elev24L = 0
    integer :: elev25L = 0
    integer :: elev26L = 0
    integer :: elev27L = 0
    integer :: elev28L = 0
    integer :: elev29L = 0
    integer :: elev30L = 0
    integer :: elev31L = 0
    integer :: elev32L = 0
    integer :: elev33L = 0
    integer :: elev34L = 0

    integer :: elev2U = 0
    integer :: elev3U = 0
    integer :: elev4U = 0
    integer :: elev5U = 0
    integer :: elev6U = 0
    integer :: elev7U = 0
    integer :: elev8U = 0
    integer :: elev9U = 0
    integer :: elev10U = 0
    integer :: elev11U = 0
    integer :: elev12U = 0
    integer :: elev13U = 0
    integer :: elev14U = 0
    integer :: elev15U = 0
    integer :: elev16U = 0
    integer :: elev17U = 0
    integer :: elev18U = 0
    integer :: elev19U = 0
    integer :: elev20U = 0
    integer :: elev23U = 0
    integer :: elev24U = 0
    integer :: elev25U = 0
    integer :: elev27U = 0
    integer :: elev28U = 0
    integer :: elev29U = 0
    integer :: elev30U = 0
    integer :: elev31U = 0
    integer :: elev33U = 0

    contains

    subroutine CreateElevationOffsets (chunk, filedatabase, qtyTemplates)
       use CFM_QuantityTemplate_m, only: CreateQtyTemplate
       use QuantityTemplates, only: AddQuantityTemplateToDatabase, &
                                    QuantityTemplate_T
       use INIT_TABLES_MODULE, only: l_elevOffset
       use Chunks_m, only: MLSChunk_T
       use MLSCommon, only: MLSFile_T

       type(MLSChunk_T), intent(in) :: chunk
       type (MLSFile_T), dimension(:), pointer :: filedatabase
       type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates

       type(QuantityTemplate_T) :: elevOffset

       ! Executables

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1A:118.B1LF:PT", &
          qName='elev1L')
       elev1L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B2LF:H2O", &
          qName='elev2L')
       elev2L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B3LF:N2O", &
          qName='elev3L')
       elev3L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B4LF:HNO3", &
          qName='elev4L')
       elev4L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B5LF:ClO", &
          qName='elev5L')
       elev5L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B6LF:O3", &
          qName='elev6L')
       elev6L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B7LF:O3", &
          qName='elev7L')
       elev7L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B8LF:PT", &
          qName='elev8L')
       elev8L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B9LF:CO", &
          qName='elev9L')
       elev9L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B10LF:ClO", &
          qName='elev10L')
       elev10L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B11LF:BrO", &
          qName='elev11L')
       elev11L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B12LF:N2O", &
          qName='elev12L')
       elev12L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B13LF:HCl", &
          qname='elev13L')
       elev13L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B14LF:O3", &
          qName='elev14L')
       elev14L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B15LF:OH", &
          qName='elev15L')
       elev15L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B16LF:OH", &
          qName='elev16L')
       elev16L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B17LF:PT", &
          qName='elev17L')
       elev17L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B18LF:OH", &
          qName='elev18L')
       elev18L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B19LF:OH", &
          qName='elev19L')
       elev19L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B20LF:PT", &
          qName='elev20L')
       elev20L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1B:118.B21LF:PT", &
          qName='elev21L')
       elev21L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1A:118.B22LD:PT", &
          qName='elev22L')
       elev22L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B23LD:H2O", &
          qName='elev23L')
       elev23L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B24LD:O3", &
          qName='elev24L')
       elev24L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B25LD:CO", &
          qName='elev25L')
       elev25L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1B:118.B26LD:PT", &
          qName='elev26L')
       elev26L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B27LM:HCN", &
          qName='elev27L')
       elev27L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B28LM:HO2", &
          qName='elev28L')
       elev28L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B29LM:HOCl", &
          qName='elev29L')
       elev29L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B30LM:HO2", &
          qName='elev30L')
       elev30L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B31LM:BrO", &
          qName='elev31L')
       elev31L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1A:118.B32LW:PT", &
          qName='elev32L')
       elev32L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B33LW:O3", &
          qName='elev33L')
       elev33L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R1B:118.B34LW:PT", &
          qName='elev34L')
       elev34L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       ! Now the the upper band
       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B2UF:H2O", &
          qName='elev2U')
       elev2U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B3UF:N2O", &
          qName='elev3U')
       elev3U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B4UF:HNO3", &
          qname='elev4U')
       elev4U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B5UF:ClO", &
          qName='elev5U')
       elev5U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B6UF:O3", &
          qName='elev6U')
       elev6U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B7UF:O3", &
          qName='elev7U')
       elev7U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B8UF:PT", &
          qName='elev8U')
       elev8U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B9UF:CO", &
          qName='elev9U')
       elev9U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B10UF:ClO", &
          qName='elev10U')
       elev10U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B11UF:BrO", &
          qName='elev11U')
       elev11U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B12UF:N2O", &
          qName='elev12U')
       elev12U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B13UF:HCl", &
          qName='elev13U')
       elev13U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B14UF:O3", &
          qname='elev14U')
       elev14U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B15UF:OH", &
          qname='elev15U')
       elev15U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B16UF:OH", &
          qName='elev16U')
       elev16U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5H:2T5.B17UF:PT", &
          qName='elev17U')
       elev17U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B18UF:OH", &
          qName='elev18U')
       elev18U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B19UF:OH", &
          qName='elev19U')
       elev19U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R5V:2T5.B20UF:PT", &
          qName='elev20U')
       elev20U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B23UD:H2O", &
          qname='elev23U')
       elev23U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B24UD:O3", &
          qName='elev24U')
       elev24U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B25UD:CO", &
          qName='elev25U')
       elev25U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R2:190.B27UM:HCN", &
          qName='elev27U')
       elev27U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B28UM:HO2", &
          qName='elev28U')
       elev28U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B29UM:HOCl", &
          qName='elev29U')
       elev29U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B30UM:HO2", &
          qName='elev30U')
       elev30U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R4:640.B31UM:BrO", &
          qName='elev31U')
       elev31U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal="R3:240.B33UW:O3", &
          qName='elev33U')
       elev33U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

    end subroutine

    subroutine FillElevationOffsets (stateVectorExtra)
       use VectorsModule, only: Vector_T, VectorValue_T, &
                                GetVectorQtyByTemplateIndex
       use CFM_Fill_M, only: SpreadFillVectorQuantity, ExplicitFillVectorQuantity

       type (Vector_T), intent(in) :: stateVectorExtra
       type(VectorValue_T) :: quantity

       ! Executables
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev1L)
       call SpreadFillVectorQuantity (quantity, 0.0187_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev2L)
       call SpreadFillVectorQuantity (quantity, 0.0151_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev3L)
       call SpreadFillVectorQuantity (quantity, 0.0150_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev4L)
       call SpreadFillVectorQuantity (quantity, 0.0148_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev5L)
       call SpreadFillVectorQuantity (quantity, 0.0146_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev6L)
       call SpreadFillVectorQuantity (quantity, 0.0144_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev7L)
       call SpreadFillVectorQuantity (quantity, -0.0002_r8)
       !print *, "Elevation offset 7L"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev8L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev9L)
       call SpreadFillVectorQuantity (quantity, -0.0004_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev10L)
       call SpreadFillVectorQuantity (quantity, 0.0077_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev11L)
       call SpreadFillVectorQuantity (quantity, 0.0077_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev12L)
       call SpreadFillVectorQuantity (quantity, 0.0077_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev13L)
       call SpreadFillVectorQuantity (quantity, 0.0076_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev14L)
       call SpreadFillVectorQuantity (quantity, 0.0076_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev15L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev16L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev17L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev18L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev19L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev20L)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev21L)
       call SpreadFillVectorQuantity (quantity, 0.0056_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev22L)
       call SpreadFillVectorQuantity (quantity, 0.0187_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev23L)
       call SpreadFillVectorQuantity (quantity, 0.0151_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev24L)
       call SpreadFillVectorQuantity (quantity, -0.0002_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev25L)
       call SpreadFillVectorQuantity (quantity, -0.0004_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev26L)
       call SpreadFillVectorQuantity (quantity, 0.0056_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev27L)
       call SpreadFillVectorQuantity (quantity, 0.0144_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev28L)
       call SpreadFillVectorQuantity (quantity, 0.0077_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev29L)
       call SpreadFillVectorQuantity (quantity, 0.0077_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev30L)
       call SpreadFillVectorQuantity (quantity, 0.0076_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev31L)
       call SpreadFillVectorQuantity (quantity, 0.0076_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev32L)
       call ExplicitFillVectorQuantity (quantity, velev32L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev33L)
       call ExplicitFillVectorQuantity (quantity, velev33L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev34L)
       call ExplicitFillVectorQuantity (quantity, velev34L)

       ! The upper bands
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev2U)
       call SpreadFillVectorQuantity (quantity, 0.0141_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev3U)
       call SpreadFillVectorQuantity (quantity, 0.0140_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev4U)
       call SpreadFillVectorQuantity (quantity, 0.0138_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev5U)
       call SpreadFillVectorQuantity (quantity, 0.0133_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev6U)
       call SpreadFillVectorQuantity (quantity, 0.0130_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev7U)
       call SpreadFillVectorQuantity (quantity, 0.0001_r8)
       !print *, "Elevation offset 7U"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev8U)
       call SpreadFillVectorQuantity (quantity, 0.0003_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev9U)
       call SpreadFillVectorQuantity (quantity, -0.0003_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev10U)
       call SpreadFillVectorQuantity (quantity, 0.0078_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev11U)
       call SpreadFillVectorQuantity (quantity, 0.0081_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev12U)
       call SpreadFillVectorQuantity (quantity, 0.0089_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev13U)
       call SpreadFillVectorQuantity (quantity, 0.0079_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev14U)
       call SpreadFillVectorQuantity (quantity, 0.0078_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev15U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev16U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev17U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev18U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev19U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev20U)
       call SpreadFillVectorQuantity (quantity, -0.0000_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev23U)
       call SpreadFillVectorQuantity (quantity, 0.0141_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev24U)
       call SpreadFillVectorQuantity (quantity, 0.0001_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev25U)
       call SpreadFillVectorQuantity (quantity, -0.0003_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev27U)
       call SpreadFillVectorQuantity (quantity, 0.0130_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev28U)
       call SpreadFillVectorQuantity (quantity, 0.0079_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev29U)
       call SpreadFillVectorQuantity (quantity, 0.0080_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev30U)
       call SpreadFillVectorQuantity (quantity, 0.0078_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev31U)
       call SpreadFillVectorQuantity (quantity, 0.0078_r8)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev33U)
       call ExplicitFillVectorQuantity (quantity, velev33U)
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
! Revision 1.5  2010/06/29 17:02:47  honghanh
! Change the identifier 'fakeChunk' to 'chunk' because
! since it is created with ChunkDivide, it's as real as a chunk
! can get.
!
! Revision 1.4  2010/06/29 16:40:23  honghanh
! Remove all function/subroutine and user type forwarding from
! all CFM modules except for from cfm.f90
!
! Revision 1.3  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
