! Copyright 2010, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
module CFM_FWDMDL_M
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
    use ForwardModelConfig, only: ForwardModelConfig_T
    use FillUtils_1, only: AutoFillVector
    use Init_tables_module, only: l_isotoperatio
    use CFM_QuantityTemplate_m, only: CreateQtyTemplate
    use QuantityTemplates, only: QuantityTemplate_T, &
                                 DestroyQuantityTemplateContents
    use VectorsModule, only: VectorValue_T, Vector_T
    use CFM_Vector_m, only: CreateValue4AgileVector, CloneAgileVector, &
                            AddValue2Vector, DestroyAgileVectorContent
    use MatrixModule_1, only: MATRIX_T
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    implicit none

    public :: ForwardModel, ForwardModel2

    interface ForwardModel
        module procedure ForwardModelObsolete
    end interface

    interface ForwardModel2
        module procedure ForwardModelWChunk, ForwardModelWMaf
    end interface

    private 
!---------------------------- RCS Ident Info -------------------------------
    character(len=*), parameter :: ModuleName= &
        "$RCSfile$"
    private :: not_used_here
!---------------------------------------------------------------------------

    type(QuantityTemplate_T), dimension(:), pointer :: qtemplates => NULL()

    contains

    ! Compute radiances based on atmospheric state information provided
    ! by fwdModelIn and fwdModelExtra
    ! This subroutine should be obsolete. Please use ForwardModelWChunk,
    ! and ForwardModelWMaf.
    subroutine ForwardModelObsolete (chunk, Config, FwdModelIn, FwdModelExtra, &
                                     FwdModelOut, Jacobian, RequestedMAF )
        use Chunks_m, only: MLSChunk_T
        use ForwardModelWrappers, only: ForwardModelOrig => ForwardModel

        ! the chunk carries the MAF to compute over
        type(MLSChunk_T), intent(in) :: chunk
        ! Configuration information for the forward model
        type(ForwardModelConfig_T), dimension(:), intent(inout) :: CONFIG
        ! Atmospheric input
        type(vector_T), intent(inout) ::  FWDMODELIN, FwdModelExtra
        ! This is the output vector where radiances are stored
        type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
        ! The matrix equivalent of the forward model with the specific
        ! setting in the config object.
        type(matrix_T), intent(inout), optional :: JACOBIAN
        ! The 0-based index of the instance (maf or profile) stored in
        ! fwdModelIn to run the forward model over.
        ! If there is only 1 instance, then this field should be 0.
        ! If there is N instance, then the appropriate values for this
        ! field is 0 to N-1.
        integer, intent(in), optional :: REQUESTEDMAF

        type(forwardModelStatus_t) :: FMSTAT ! Reverse comm. stuff
        integer :: i
        integer :: FinalMAF

        fmStat%newScanHydros = .true.
        if (present(jacobian)) then
            call allocate_test(fmstat%rows, jacobian%row%nb, "fmStat%rows", moduleName)
        end if
      
        do i=1, size(config)
            if ( present ( requestedMAF ) ) then
                fmStat%maf = requestedMAF
                finalMAF = requestedMAF
            else
                fmStat%maf = 0
                finalMAF = chunk%lastMAFIndex - chunk%firstMAFIndex
            end if
            ! Loop over MAFs
            do while (fmStat%maf <= finalMAF )
                fmStat%maf = fmStat%maf + 1
                call ForwardModelOrig (config(i), fwdmodelIn, fwdModelExtra, fwdModelOut, &
                fmStat, Jacobian)
            end do
        end do

        if (present(jacobian)) then
            call deallocate_test(fmStat%rows, "fmStat%rows", moduleName)
        end if

    end subroutine

    subroutine PreForwardModel (config, fwdModelExtra, fmStat, jacobian)
        type(ForwardModelConfig_T), dimension(:), intent(in) :: CONFIG
        type(vector_T), intent(inout) ::  FwdModelExtra
        type(forwardModelStatus_t), intent(inout) :: FMSTAT
        type(matrix_T), intent(in), optional :: JACOBIAN

        type(QuantityTemplate_T) :: qtemplate
        type(VectorValue_T) :: value
        integer :: i, b, m, sideband
        integer :: counter

        counter = 0
        ! Insert isotope quantities into extra vector
        sideband = 1
        do i = 1, size(config)
            do b = 1, size(config(i)%beta_group)
                if (config(i)%beta_group(b)%group) then ! A molecule group
                    ! First LBL molecules' ratios
                    do m = 1, size(config(i)%beta_group(b)%lbl(sideband)%molecules)
                        qtemplate = CreateQtyTemplate(l_isotoperatio, &
                        qmolecule=config(i)%beta_group(b)%lbl(sideband)%molecules(m))
                        value = CreateValue4AgileVector(qtemplate)
                        call AddValue2Vector(fwdModelExtra, value)
                        counter = counter + 1
                    end do !m
                    if (associated(config(i)%beta_group(b)%pfa(sideband)%molecules)) then
                        ! Now PFA molecules' ratios
                        do m = 1, size(config(i)%beta_group(b)%pfa(sideband)%molecules)
                            qtemplate = CreateQtyTemplate(l_isotoperatio, &
                            qmolecule=config(i)%beta_group(b)%pfa(sideband)%molecules(m))
                            value = CreateValue4AgileVector(qtemplate)
                            call AddValue2Vector(fwdModelExtra, value)
                            counter = counter + 1
                        end do !m
                    end if
                end if
            end do ! b
        end do ! i
        ! Fill isotope in the fwdModelExtra vector
        call AutoFillVector(fwdModelExtra)
        ! Record newly created template object
        allocate(qtemplates(counter), stat=b)
        if (b /= 0) &
            call MLSMessage (MLSMSG_Error, moduleName, &
            "Cannot allocate qtemplates")
        m = size(fwdModelExtra%quantities) - counter + 1
        do i = 1, counter
            qtemplates(i) = fwdModelExtra%quantities(m)%template
            m = m + 1
        end do
 
        fmStat%newScanHydros = .true.
        if (present(jacobian)) then
            call allocate_test(fmstat%rows, jacobian%row%nb, "fmStat%rows", moduleName)
        else
            call allocate_test(fmstat%rows, 0, "fmStat%rows", moduleName)
        end if
    end subroutine

    subroutine PostForwardModel (fmStat)
        type(forwardModelStatus_t), intent(inout) :: FMSTAT
        integer :: i

        do i = 1, size(qtemplates)
            call DestroyQuantityTemplateContents(qtemplates(i))
        end do
        deallocate(qtemplates)
        call deallocate_test(fmStat%rows, "fmStat%rows", moduleName)
    end subroutine

    ! This method is obsolete.
    ! Compute radiances based on atmospheric state information provided
    ! by fwdModelIn and fwdModelExtra.
    subroutine ForwardModelWChunk (chunk, Config, FwdModelIn, FwdModelExtra, &
                                    FwdModelOut, Jacobian )
        use Chunks_m, only: MLSChunk_T
        use ForwardModelWrappers, only: ForwardModelOrig => ForwardModel

        ! the chunk carries the MAF to compute over
        type(MLSChunk_T), intent(in) :: chunk
        ! Configuration information for the forward model
        type(ForwardModelConfig_T), dimension(:), intent(inout) :: CONFIG
        ! Atmospheric input
        type(vector_T), intent(inout) ::  FWDMODELIN, FwdModelExtra
        ! This is the output vector where radiances are stored
        type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
        ! The matrix equivalent of the forward model with the specific
        ! setting in the config object.
        type(matrix_T), intent(inout), optional :: JACOBIAN

        type(forwardModelStatus_t) :: FMSTAT ! Reverse comm. stuff
        type(Vector_T) :: adjustedExtraInput
        integer :: FinalMAF
        integer :: i

        call cloneAgileVector(adjustedExtraInput, fwdModelExtra)
        call preForwardModel (config, adjustedExtraInput, fmStat, jacobian)

        do i=1, size(config)
            fmStat%maf = 0
            finalMAF = chunk%lastMAFIndex - chunk%firstMAFIndex

            ! Loop over MAFs
            do while (fmStat%maf <= finalMAF )
                fmStat%maf = fmStat%maf + 1
                call ForwardModelOrig (config(i), fwdmodelIn, &
                adjustedExtraInput, fwdModelOut, fmStat, Jacobian)
            end do
        end do

        call DestroyAgileVectorContent(adjustedExtraInput)
        call postForwardModel(fmStat)
    end subroutine

    ! Compute radiances based on atmospheric state information provided
    ! by fwdModelIn and fwdModelExtra.
    subroutine ForwardModelWMaf (RequestedMAF, Config, FwdModelIn, FwdModelExtra, &
                                    FwdModelOut, Jacobian )
        use ForwardModelWrappers, only: ForwardModelOrig => ForwardModel
        use ForwardModelVectorTools
        use init_tables_module, only: l_ptan

        ! The 0-based index of the instance (maf or profile) stored in
        ! fwdModelIn to run the forward model over.
        ! If there is only 1 instance, then this field should be 0.
        ! If there is N instance, then the appropriate values for this
        ! field is 0 to N-1.
        integer, intent(in) :: RequestedMAF
        ! Configuration information for the forward model
        type(ForwardModelConfig_T), dimension(:), intent(inout) :: CONFIG
        ! Atmospheric input
        type(vector_T), intent(inout) ::  FWDMODELIN, FwdModelExtra
        ! This is the output vector where radiances are stored
        type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
        ! The matrix equivalent of the forward model with the specific
        ! setting in the config object.
        type(matrix_T), intent(inout), optional :: JACOBIAN

        type(forwardModelStatus_t) :: FMSTAT ! Reverse comm. stuff
        type(Vector_T) :: adjustedExtraInput
        integer :: i
        type(VectorValue_T), pointer :: ptan

        call cloneagilevector(adjustedExtraInput, fwdModelExtra)
        call preForwardModel (config, adjustedExtraInput, fmStat, jacobian)

        fmStat%maf = requestedMaf + 1
        do i=1, size(config)
            call ForwardModelOrig (config(i), fwdmodelIn, &
            adjustedExtraInput, fwdModelOut, fmStat, Jacobian)
        end do

        call DestroyAgileVectorContent(adjustedExtraInput)
        call postForwardModel(fmStat)
        
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
! Revision 1.11  2011/12/24 18:39:13  honghanh
! Clean up unused imports and variables
!
! Revision 1.10  2011/12/23 22:56:12  honghanh
! Add AutoFillVector call to ForwardModel2 subroutine,
! to automatically add and fill isotope ratio in beta group
! (according to L2CF configuration file)
!
! Revision 1.9  2011/12/15 18:27:44  honghanh
! Documentation and code clean up, including removing unused and broken
! subroutines.
!
! Revision 1.8  2011/06/21 18:54:28  honghanh
! Obsoleting the old ForwardModel subroutine.
! Instead, adding 2 new subroutines without using requested maf as an optional parameter.
! Introducing ForwardModel2
!
! Revision 1.7  2010/11/19 17:16:08  honghanh
! Add example to use the requestedMAF feature in forward model
!
! Revision 1.6  2010/11/19 01:04:51  livesey
! Added the optional RequestdMAF argument
!
! Revision 1.5  2010/08/05 16:23:03  honghanh
! Added Jacobian to forwardModel subroutine
!
! Revision 1.4  2010/06/29 17:02:47  honghanh
! Change the identifier 'fakeChunk' to 'chunk' because
! since it is created with ChunkDivide, it's as real as a chunk
! can get.
!
! Revision 1.3  2010/06/29 15:53:45  honghanh
! Add copyright comments and support for CVS log in the file
!
