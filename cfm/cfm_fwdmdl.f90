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

   implicit none

   public :: ForwardModel

!---------------------------- RCS Ident Info -------------------------------
   character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
   private :: not_used_here
!---------------------------------------------------------------------------

   private

   contains

   ! Compute radiances based on atmospheric state information provided
   ! by fwdModelIn and fwdModelExtra
   subroutine ForwardModel (chunk, Config, FwdModelIn, FwdModelExtra, &
                            FwdModelOut, Jacobian, RequestedMAF )
      use VectorsModule, only: Vector_T
      use MatrixModule_1, only: MATRIX_T
      use Chunks_m, only: MLSChunk_T
      use ForwardModelWrappers, only: ForwardModelOrig => ForwardModel
      use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

      ! the chunk carries the MAF to compute over
      type(MLSChunk_T), intent(in) :: chunk
      ! Configuration information for the forward model
      type(ForwardModelConfig_T), dimension(:), intent(inout) :: CONFIG
      ! Atmospheric input
      type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
      ! This is the output vector where radiances are stored
      type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
      ! The matrix equivalent of the forward model with the specific
      ! setting in the config object.
      type(matrix_T), intent(inout), optional :: JACOBIAN
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
