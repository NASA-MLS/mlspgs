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
    use CFM_Constants_m

    implicit none

    private

    public :: CreateElevationOffsets, FillElevationOffsets

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

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
          filedatabase=filedatabase, chunk=chunk, qSignal=band1L, &
          qName='elev1L')
       elev1L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band2L, &
          qName='elev2L')
       elev2L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band3L, &
          qName='elev3L')
       elev3L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band4L, &
          qName='elev4L')
       elev4L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band5L, &
          qName='elev5L')
       elev5L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band6L, &
          qName='elev6L')
       elev6L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band7L, &
          qName='elev7L')
       elev7L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band8L, &
          qName='elev8L')
       elev8L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band9L, &
          qName='elev9L')
       elev9L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band10L, &
          qName='elev10L')
       elev10L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band11L, &
          qName='elev11L')
       elev11L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band12L, &
          qName='elev12L')
       elev12L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band13L, &
          qname='elev13L')
       elev13L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band14L, &
          qName='elev14L')
       elev14L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band15L, &
          qName='elev15L')
       elev15L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band16L, &
          qName='elev16L')
       elev16L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band17L, &
          qName='elev17L')
       elev17L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band18L, &
          qName='elev18L')
       elev18L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band19L, &
          qName='elev19L')
       elev19L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band20L, &
          qName='elev20L')
       elev20L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band21L, &
          qName='elev21L')
       elev21L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band22L, &
          qName='elev22L')
       elev22L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band23L, &
          qName='elev23L')
       elev23L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band24L, &
          qName='elev24L')
       elev24L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band25L, &
          qName='elev25L')
       elev25L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band26L, &
          qName='elev26L')
       elev26L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band27L, &
          qName='elev27L')
       elev27L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band28L, &
          qName='elev28L')
       elev28L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band29L, &
          qName='elev29L')
       elev29L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band30L, &
          qName='elev30L')
       elev30L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band31L, &
          qName='elev31L')
       elev31L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band32L, &
          qName='elev32L')
       elev32L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band33L, &
          qName='elev33L')
       elev33L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band34L, &
          qName='elev34L')
       elev34L = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       ! Now the the upper band
       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band2U, &
          qName='elev2U')
       elev2U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band3U, &
          qName='elev3U')
       elev3U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band4U, &
          qname='elev4U')
       elev4U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band5U, &
          qName='elev5U')
       elev5U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band6U, &
          qName='elev6U')
       elev6U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band7U, &
          qName='elev7U')
       elev7U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band8U, &
          qName='elev8U')
       elev8U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band9U, &
          qName='elev9U')
       elev9U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band10U, &
          qName='elev10U')
       elev10U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band11U, &
          qName='elev11U')
       elev11U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band12U, &
          qName='elev12U')
       elev12U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band13U, &
          qName='elev13U')
       elev13U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band14U, &
          qname='elev14U')
       elev14U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band15U, &
          qname='elev15U')
       elev15U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band16U, &
          qName='elev16U')
       elev16U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band17U, &
          qName='elev17U')
       elev17U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band18U, &
          qName='elev18U')
       elev18U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band19U, &
          qName='elev19U')
       elev19U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band20U, &
          qName='elev20U')
       elev20U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band23U, &
          qname='elev23U')
       elev23U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band24U, &
          qName='elev24U')
       elev24U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band25U, &
          qName='elev25U')
       elev25U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band27U, &
          qName='elev27U')
       elev27U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band28U, &
          qName='elev28U')
       elev28U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band29U, &
          qName='elev29U')
       elev29U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band30U, &
          qName='elev30U')
       elev30U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band31U, &
          qName='elev31U')
       elev31U = AddQuantityTemplateToDatabase(qtyTemplates, elevOffset)

       elevOffset = CreateQtyTemplate(l_elevOffset, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band33U, &
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
       call SpreadFillVectorQuantity (quantity, velev1L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev2L)
       call SpreadFillVectorQuantity (quantity, velev2L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev3L)
       call SpreadFillVectorQuantity (quantity, velev3L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev4L)
       call SpreadFillVectorQuantity (quantity, velev4L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev5L)
       call SpreadFillVectorQuantity (quantity, velev5L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev6L)
       call SpreadFillVectorQuantity (quantity, velev6L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev7L)
       call SpreadFillVectorQuantity (quantity, velev7L)
       !print *, "Elevation offset 7L"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev8L)
       call SpreadFillVectorQuantity (quantity, velev8L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev9L)
       call SpreadFillVectorQuantity (quantity, velev9L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev10L)
       call SpreadFillVectorQuantity (quantity, velev10L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev11L)
       call SpreadFillVectorQuantity (quantity, velev11L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev12L)
       call SpreadFillVectorQuantity (quantity, velev12L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev13L)
       call SpreadFillVectorQuantity (quantity, velev13L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev14L)
       call SpreadFillVectorQuantity (quantity, velev14L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev15L)
       call SpreadFillVectorQuantity (quantity, velev15L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev16L)
       call SpreadFillVectorQuantity (quantity, velev16L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev17L)
       call SpreadFillVectorQuantity (quantity, velev17L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev18L)
       call SpreadFillVectorQuantity (quantity, velev18L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev19L)
       call SpreadFillVectorQuantity (quantity, velev19L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev20L)
       call SpreadFillVectorQuantity (quantity, velev20L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev21L)
       call SpreadFillVectorQuantity (quantity, velev21L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev22L)
       call SpreadFillVectorQuantity (quantity, velev22L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev23L)
       call SpreadFillVectorQuantity (quantity, velev23L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev24L)
       call SpreadFillVectorQuantity (quantity, velev24L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev25L)
       call SpreadFillVectorQuantity (quantity, velev25L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev26L)
       call SpreadFillVectorQuantity (quantity, velev26L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev27L)
       call SpreadFillVectorQuantity (quantity, velev27L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev28L)
       call SpreadFillVectorQuantity (quantity, velev28L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev29L)
       call SpreadFillVectorQuantity (quantity, velev29L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev30L)
       call SpreadFillVectorQuantity (quantity, velev30L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev31L)
       call SpreadFillVectorQuantity (quantity, velev31L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev32L)
       call ExplicitFillVectorQuantity (quantity, velev32L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev33L)
       call ExplicitFillVectorQuantity (quantity, velev33L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev34L)
       call ExplicitFillVectorQuantity (quantity, velev34L)

       ! The upper bands
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev2U)
       call SpreadFillVectorQuantity (quantity, velev2U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev3U)
       call SpreadFillVectorQuantity (quantity, velev3U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev4U)
       call SpreadFillVectorQuantity (quantity, velev4U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev5U)
       call SpreadFillVectorQuantity (quantity, velev5U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev6U)
       call SpreadFillVectorQuantity (quantity, velev6U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev7U)
       call SpreadFillVectorQuantity (quantity, velev7U)
       !print *, "Elevation offset 7U"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev8U)
       call SpreadFillVectorQuantity (quantity, velev8U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev9U)
       call SpreadFillVectorQuantity (quantity, velev9U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev10U)
       call SpreadFillVectorQuantity (quantity, velev10U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev11U)
       call SpreadFillVectorQuantity (quantity, velev11U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev12U)
       call SpreadFillVectorQuantity (quantity, velev12U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev13U)
       call SpreadFillVectorQuantity (quantity, velev13U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev14U)
       call SpreadFillVectorQuantity (quantity, velev14U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev15U)
       call SpreadFillVectorQuantity (quantity, velev15U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev16U)
       call SpreadFillVectorQuantity (quantity, velev16U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev17U)
       call SpreadFillVectorQuantity (quantity, velev17U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev18U)
       call SpreadFillVectorQuantity (quantity, velev18U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev19U)
       call SpreadFillVectorQuantity (quantity, velev19U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev20U)
       call SpreadFillVectorQuantity (quantity, velev20U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev23U)
       call SpreadFillVectorQuantity (quantity, velev23U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev24U)
       call SpreadFillVectorQuantity (quantity, velev24U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev25U)
       call SpreadFillVectorQuantity (quantity, velev25U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev27U)
       call SpreadFillVectorQuantity (quantity, velev27U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev28U)
       call SpreadFillVectorQuantity (quantity, velev28U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev29U)
       call SpreadFillVectorQuantity (quantity, velev29U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev30U)
       call SpreadFillVectorQuantity (quantity, velev30U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, elev31U)
       call SpreadFillVectorQuantity (quantity, velev31U)
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
! Revision 1.6  2010/11/03 18:34:46  honghanh
! Add qName as an optional argument to CreateQtyTemplate.
! This is to help debugging process.
!
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
