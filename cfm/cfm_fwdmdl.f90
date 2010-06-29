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
   subroutine ForwardModel (fakeChunk, Config, FwdModelIn, FwdModelExtra, &
                            FwdModelOut, Jacobian)
      use VectorsModule, only: Vector_T
      use MatrixModule_1, only: MATRIX_T
      use Chunks_m, only: MLSChunk_T
      use ForwardModelWrappers, only: ForwardModelOrig => ForwardModel

      ! the chunk carries the MAF to compute over
      type(MLSChunk_T), intent(in) :: fakeChunk
      ! Configuration information for the forward model
      type(ForwardModelConfig_T), dimension(:), intent(inout) :: CONFIG
      ! Atmospheric input
      type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
      ! This is the output vector where radiances are stored
      type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
      ! The matrix equivalent of the forward model with the specific
      ! setting in the config object.
      type(matrix_T), intent(inout), optional :: JACOBIAN

      type(forwardModelStatus_t) :: FMSTAT ! Reverse comm. stuff
      integer :: i

      fmStat%newScanHydros = .true.

      do i=1, size(config)
         fmStat%maf = 0
         ! Loop over MAFs
         do while (fmStat%maf < fakeChunk%lastMAFIndex-fakeChunk%firstMAFIndex+1)
            fmStat%maf = fmStat%maf + 1
            call ForwardModelOrig (config(i), fwdmodelIn, fwdModelExtra, fwdModelOut, &
                                   fmStat, Jacobian)
         end do
      end do

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
