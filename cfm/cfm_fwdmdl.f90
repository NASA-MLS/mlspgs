module CFM_FWDMDL_M
   use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
   use ForwardModelConfig, only: ForwardModelConfig_T

   implicit none

   public :: ForwardModel
   public :: ForwardModelConfig_T, FORWARDMODELSTATUS_T

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
end module
