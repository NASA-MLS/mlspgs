! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module cfm_lsf_m ! This module is to help cfm_mlssetup subroutine.
                 ! Quantities allocated in here should be deallocate by cfm_mlssetup_m
    use MLSCommon, only: r8
    use CFM_Constants_m

    implicit none

    private

    public :: CreateLimbSidebandFractions, FillLimbSidebandFractions

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

    integer :: limbSidebandFraction1L = 0
    integer :: limbSidebandFraction2L = 0
    integer :: limbSidebandFraction3L = 0
    integer :: limbSidebandFraction4L = 0
    integer :: limbSidebandFraction5L = 0
    integer :: limbSidebandFraction6L = 0
    integer :: limbSidebandFraction7L = 0
    integer :: limbSidebandFraction8L = 0
    integer :: limbSidebandFraction9L = 0
    integer :: limbSidebandFraction10L = 0
    integer :: limbSidebandFraction11L = 0
    integer :: limbSidebandFraction12L = 0
    integer :: limbSidebandFraction13L = 0
    integer :: limbSidebandFraction14L = 0
    integer :: limbSidebandFraction15L = 0
    integer :: limbSidebandFraction16L = 0
    integer :: limbSidebandFraction17L = 0
    integer :: limbSidebandFraction18L = 0
    integer :: limbSidebandFraction19L = 0
    integer :: limbSidebandFraction20L = 0
    integer :: limbSidebandFraction21L = 0
    integer :: limbSidebandFraction22L = 0
    integer :: limbSidebandFraction23L = 0
    integer :: limbSidebandFraction24L = 0
    integer :: limbSidebandFraction25L = 0
    integer :: limbSidebandFraction26L = 0
    integer :: limbSidebandFraction27L = 0
    integer :: limbSidebandFraction28L = 0
    integer :: limbSidebandFraction29L = 0
    integer :: limbSidebandFraction30L = 0
    integer :: limbSidebandFraction31L = 0
    integer :: limbSidebandFraction32L = 0
    integer :: limbSidebandFraction33L = 0
    integer :: limbSidebandFraction34L = 0

    integer :: limbSidebandFraction2U = 0
    integer :: limbSidebandFraction3U = 0
    integer :: limbSidebandFraction4U = 0
    integer :: limbSidebandFraction5U = 0
    integer :: limbSidebandFraction6U = 0
    integer :: limbSidebandFraction7U = 0
    integer :: limbSidebandFraction8U = 0
    integer :: limbSidebandFraction9U = 0
    integer :: limbSidebandFraction10U = 0
    integer :: limbSidebandFraction11U = 0
    integer :: limbSidebandFraction12U = 0
    integer :: limbSidebandFraction13U = 0
    integer :: limbSidebandFraction14U = 0
    integer :: limbSidebandFraction15U = 0
    integer :: limbSidebandFraction16U = 0
    integer :: limbSidebandFraction17U = 0
    integer :: limbSidebandFraction18U = 0
    integer :: limbSidebandFraction19U = 0
    integer :: limbSidebandFraction20U = 0
    integer :: limbSidebandFraction23U = 0
    integer :: limbSidebandFraction24U = 0
    integer :: limbSidebandFraction25U = 0
    integer :: limbSidebandFraction27U = 0
    integer :: limbSidebandFraction28U = 0
    integer :: limbSidebandFraction29U = 0
    integer :: limbSidebandFraction30U = 0
    integer :: limbSidebandFraction31U = 0
    integer :: limbSidebandFraction33U = 0

    contains

    ! InitQuantityTemplates should be called priori to this subroutine
    subroutine CreateLimbSidebandFractions (chunk, filedatabase, qtyTemplates)
       use CFM_QuantityTemplate_m, only: CreateQtyTemplate
       use QuantityTemplates, only: AddQuantityTemplateToDatabase, &
                                    QuantityTemplate_T
       use INIT_TABLES_MODULE, only: l_limbsidebandFraction
       use Chunks_m, only: MLSChunk_T
       use MLSCommon, only: MLSFile_T

       type(MLSChunk_T), intent(in) :: chunk
       type (MLSFile_T), dimension(:), pointer :: filedatabase
       type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates

       type(QuantityTemplate_T) :: lsbFraction

       ! Executables

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band1L, &
          qName='limbSidebandFraction1L')
       limbSidebandFraction1L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band2L, &
          qName='limbSidebandFraction2L')
       limbSidebandFraction2L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band3L, &
          qName='limbSidebandFraction3L')
       limbSidebandFraction3L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band4L, &
          qName='limbSidebandFraction4L')
       limbSidebandFraction4L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band5L, &
          qName='limbSidebandFraction5L')
       limbSidebandFraction5L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band6L, &
          qName='limbSidebandFraction6L')
       limbSidebandFraction6L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band7L, &
          qName='limbSidebandFraction7L')
       limbSidebandFraction7L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band8L, &
          qName='limbSidebandFraction8L')
       limbSidebandFraction8L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band9L, &
          qName='limbSidebandFraction9L')
       limbSidebandFraction9L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band10L, &
          qName='limbSidebandFraction10L')
       limbSidebandFraction10L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band11L, &
          qName='limbSidebandFraction11L')
       limbSidebandFraction11L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band12L, &
          qName='limbSidebandFraction12L')
       limbSidebandFraction12L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band13L, &
          qName='limbSidebandFraction13L')
       limbSidebandFraction13L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band14L, &
          qName='limbSidebandFraction14L')
       limbSidebandFraction14L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band15L, &
          qName='limbSidebandFraction15L')
       limbSidebandFraction15L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band16L, &
          qName='limbSidebandFraction16L')
       limbSidebandFraction16L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band17L, &
          qName='limbSidebandFraction17L')
       limbSidebandFraction17L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band18L, &
          qName='limbSidebandFraction18L')
       limbSidebandFraction18L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band19L, &
          qName='limbSidebandFraction19L')
       limbSidebandFraction19L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band20L, &
          qName='limbSidebandFraction20L')
       limbSidebandFraction20L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band21L, &
          qName='limbSidebandFraction21L')
       limbSidebandFraction21L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band22L, &
          qName='limbSidebandFraction22L')
       limbSidebandFraction22L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band23L, &
          qName='limbSidebandFraction23L')
       limbSidebandFraction23L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band24L, &
          qName='limbSidebandFraction24L')
       limbSidebandFraction24L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band25L, &
          qName='limbSidebandFraction25L')
       limbSidebandFraction25L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band26L, &
          qName='limbSidebandFraction26L')
       limbSidebandFraction26L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band27L, &
          qName='limbSidebandFraction27L')
       limbSidebandFraction27L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band28L, &
          qName='limbSidebandFraction28L')
       limbSidebandFraction28L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band29L, &
          qName='limbSidebandFraction29L')
       limbSidebandFraction29L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band30L, &
          qName='limbSidebandFraction30L')
       limbSidebandFraction30L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band31L, &
          qName='limbSidebandFraction31L')
       limbSidebandFraction31L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band32L, &
          qName='limbSidebandFraction32L')
       limbSidebandFraction32L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band33L, &
          qName='limbSidebandFraction33L')
       limbSidebandFraction33L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band34L, &
          qName='limbSidebandFraction34L')
       limbSidebandFraction34L = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       ! Now the the upper band
       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band2U, &
          qName='limbSidebandFraction2U')
       limbSidebandFraction2U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band3U, &
          qName='limbSidebandFraction3U')
       limbSidebandFraction3U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band4U, &
          qName='limbSidebandFraction4U')
       limbSidebandFraction4U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band5U, &
          qName='limbSidebandFraction5U')
       limbSidebandFraction5U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band6U, &
          qName='limbSidebandFraction6U')
       limbSidebandFraction6U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band7U, &
          qName='limbSidebandFraction7U')
       limbSidebandFraction7U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band8U, &
          qName='limbSidebandFraction8U')
       limbSidebandFraction8U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band9U, &
          qName='limbSidebandFraction9U')
       limbSidebandFraction9U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band10U, &
          qName='limbSidebandFraction10U')
       limbSidebandFraction10U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band11U, &
          qName='limbSidebandFraction11U')
       limbSidebandFraction11U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band12U, &
          qName='limbSidebandFraction12U')
       limbSidebandFraction12U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band13U, &
          qName='limbSidebandFraction13U')
       limbSidebandFraction13U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band14U, &
          qName='limbSidebandFraction14U')
       limbSidebandFraction14U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band15U, &
          qName='limbSidebandFraction15U')
       limbSidebandFraction15U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band16U, &
          qName='limbSidebandFraction16U')
       limbSidebandFraction16U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band17U, &
          qName='limbSidebandFraction17U')
       limbSidebandFraction17U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band18U, &
          qName='limbSidebandFraction18U')
       limbSidebandFraction18U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band19U, &
          qName='limbSidebandFraction19U')
       limbSidebandFraction19U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band20U, &
          qName='limbSidebandFraction20U')
       limbSidebandFraction20U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band23U, &
          qName='limbSidebandFraction23U')
       limbSidebandFraction23U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band24U, &
          qName='limbSidebandFraction24U')
       limbSidebandFraction24U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band25U, &
          qName='limbSidebandFraction25U')
       limbSidebandFraction25U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band27U, &
          qName='limbSidebandFraction27U')
       limbSidebandFraction27U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band28U, &
          qName='limbSidebandFraction28U')
       limbSidebandFraction28U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band29U, &
          qName='limbSidebandFraction29U')
       limbSidebandFraction29U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band30U, &
          qName='limbSidebandFraction30U')
       limbSidebandFraction30U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band31U, &
          qName='limbSidebandFraction31U')
       limbSidebandFraction31U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

       lsbFraction = CreateQtyTemplate(l_limbsidebandFraction, &
          filedatabase=filedatabase, chunk=chunk, qSignal=band33U, &
          qName='limbSidebandFraction33U')
       limbSidebandFraction33U = AddQuantityTemplateToDatabase(qtyTemplates, lsbFraction)

    end subroutine

    subroutine FillLimbSidebandFractions (stateVectorExtra)
       use VectorsModule, only: Vector_T, VectorValue_T, Dump, &
                                GetVectorQtyByTemplateIndex
       use CFM_Fill_M, only: ExplicitFillVectorQuantity, SpreadFillVectorQuantity

       type (Vector_T), intent(in) :: stateVectorExtra
       type(VectorValue_T) :: quantity

       ! Executables
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction1L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction1L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction2L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction2L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction3L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction3L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction4L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction4L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction5L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction5L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction6L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction6L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction7L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction7L)
       !print *, "LimbSidebandFraction7L value"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction8L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction8L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction9L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction9L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction10L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction10L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction11L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction11L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction12L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction12L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction13L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction13L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction14L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction14L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction15L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction15L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction16L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction16L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction17L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction17L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction18L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction18L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction19L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction19L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction20L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction20L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction21L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction21L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction22L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction22L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction23L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction23L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction24L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction24L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction25L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction25L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction26L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction26L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction27L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction27L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction28L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction28L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction29L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction29L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction30L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction30L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction31L)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction31L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction32L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction32L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction33L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction33L)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction34L)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction34L)

       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction2U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction2U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction3U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction3U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction4U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction4U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction5U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction5U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction6U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction6U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction7U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction7U)
       !print *, "LimbSidebandFraction7U value"
       !call dump(quantity, details=3)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction8U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction8U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction9U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction9U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction10U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction10U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction11U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction11U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction12U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction12U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction13U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction13U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction14U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction14U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction15U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction15U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction16U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction16U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction17U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction17U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction18U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction18U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction19U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction19U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction20U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction20U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction23U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction23U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction24U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction24U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction25U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction25U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction27U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction27U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction28U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction28U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction29U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction29U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction30U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction30U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction31U)
       call SpreadFillVectorQuantity (quantity, vlimbSidebandFraction31U)
       quantity = GetVectorQtyByTemplateIndex (stateVectorExtra, limbSidebandFraction33U)
       call ExplicitFillVectorQuantity (quantity, vlimbSidebandFraction33U)

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
! Revision 1.8  2011/10/18 16:59:29  honghanh
! Removing declaration of constants since they were already refractored.
!
! Revision 1.7  2011/10/18 16:56:47  honghanh
! Refractoring limb sideband constants to CFM_Constants_m.
!
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
