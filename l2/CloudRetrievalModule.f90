!=============================================================================
module CloudRetrievalModule
!=============================================================================

!{This module inverts the radiative transfer equation, to solve for
! atmospheric parameters, given radiance measurements.
!
! This module and ones it calls consume most of the cycles.

  implicit NONE
  private
  public :: CloudRetrieval

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
       
contains

! ------------------------------------------  CloudRetrieval  -----
  subroutine CloudRetrieval(Method, ConfigDatabase,configIndices,fwdModelExtra,&
      & measurements,MeasurementSD, state, OutputSD, Covariance, &
      & jacobian, chunk,maxJacobians,initlambda)

      use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
                                        & MLSMSG_Warning, MLSMSG_Deallocate
      use MLSCommon, only: MLSCHUNK_T, R8, RM, RV
      use Intrinsic, only: L_PTAN, L_RADIANCE,                           &
                     & L_CLOUDINDUCEDRADIANCE,                           &
                     & L_CLOUDEXTINCTION,                                &
                     & L_CLOUDRADSENSITIVITY
      use Init_Tables_Module, only: L_highcloud, L_lowcloud, l_lostransfunc, &
         & l_clear, l_earthradius
      use MatrixModule_0, only: MatrixInversion, MATRIXELEMENT_T
      use MatrixModule_1, only: AddToMatrixDatabase, CreateEmptyMatrix, ClearMatrix,&
      & DestroyMatrix, FINDBLOCK, GetDiagonal, GetFromMatrixDatabase, Matrix_T, &
      & Matrix_Database_T, Matrix_SPD_T, MultiplyMatrixVectorNoT, &
      & Sparsify, MultiplyMatrix_XTY, UpdateDiagonal
      use MLSSignals_m, only: SIGNAL_T
      use ForwardModelConfig, only: ForwardModelConfig_T
      use ForwardModelWrappers, only: ForwardModel
      use ForwardModelIntermediate, only: ForwardModelIntermediate_T, &
        & ForwardModelStatus_T
      use VectorsModule, only: ClearMask, ClearUnderMask, &
       & ClearVector, CloneVector, CopyVector, CopyVectorMask, CreateMask, &
       & DestroyVectorInfo, DumpMask, GetVectorQuantityByType, &
       & IsVectorQtyMasked, M_Fill, M_FullDerivatives, M_LinAlg, M_Cloud, &
       & M_Tikhonov, Vector_T, VectorValue_T

      ! IO variables
      type(vector_T), pointer :: Apriori  ! A priori estimate of state
      type(forwardModelConfig_T), dimension(:), pointer :: ConfigDatabase
      integer, pointer, dimension(:) :: ConfigIndices    ! In ConfigDatabase
      type(matrix_SPD_T), pointer :: Covariance     ! covariance**(-1) of Apriori
      logical :: Diagonal                 ! "Iterate with the diagonal of the
                                       ! a priori covariance matrix until
                                       ! convergence, then put in the whole
                                       ! thing and iterate until it converges
                                       ! again (hopefully only once).
      type(vector_T), pointer :: FwdModelExtra
      type(MLSChunk_T) :: CHUNK
      real(r8) :: InitLambda              ! Initial Levenberg-Marquardt parameter
      integer :: MaxJacobians             ! Maximum number of Jacobian
                                          ! evaluations of Newtonian method
      integer :: Method                   ! Method to use for inversion
      type(vector_T), pointer :: Measurements  ! The measurements vector
      type(vector_T), pointer :: MeasurementSD ! The measurements vector's Std. Dev.
      type(vector_T), pointer :: OutputSD ! Vector containing SD of result
      type(vector_T), pointer :: State    ! The state vector
      type(matrix_T), pointer :: Jacobian ! The Jacobian matrix

      ! Local variables
      type (ForwardModelStatus_T) :: FmStat        ! Status for forward model
      type (vector_T) :: FwdModelOut1               ! Forward outputs
      type (Signal_T) :: signal                    ! signal info in each model

      type(MatrixElement_T), pointer :: JBLOCK       ! A block from the jacobian
      type(vector_T) :: CovarianceDiag  ! Diagonal of apriori Covariance  
      
      integer :: coljBlock     ! Column index for jacobian
      integer :: rowjBlock     ! Row index for jacobian
      integer :: nMAFs                    ! number of mafs
      integer :: nFreqs      ! number of frequencies in each block
      real(r8) :: badValue

      ! ------------------------------------------      
      select case (method)
         case (l_lowcloud)
            call LowCloudRetrieval
         case (l_highcloud)
            call HighCloudRetrieval
      end select


  contains      
  ! ------------------------------------------  HighCloudRetrieval  -----
    subroutine HighCloudRetrieval
      
      ! Local Variables
      type (ForwardModelIntermediate_T) :: Fmw     ! Work space for forward model
      type (VectorValue_T), pointer :: xExtPtr        ! pointer of l_cloudExtinction quantity
      type (VectorValue_T), pointer :: xExtVar        ! variance of apriori
      type (VectorValue_T), pointer :: outExtSD        ! SD of output
      type (VectorValue_T), pointer :: ConstrainTcir        ! for model cloud top indicator
      type (VectorValue_T), pointer :: Tcir        ! cloud-induced radiance
      type (VectorValue_T), pointer :: Terr        ! cloud-induced radiance SD
      type (VectorValue_T), pointer :: PTAN        ! Tgt pressure
!      type (VectorValue_T), pointer :: Re        ! Earth Radius
      type (VectorValue_T), pointer :: Tb0         ! model clear sky radiance (100%RH)
      type (VectorValue_T), pointer :: Slope       ! sensitivity slope to convert cloud
                                                   ! radiance to optical depth
                                          
      integer :: i,j,k,i1,j1,ich,maf,mif          ! Loop subscripts
      integer :: nz        ! number of retrieval levels
      integer :: nInst     ! number of retrieval instances in the chunk
      integer :: cloudysky   ! cloudysky index from Model Configuration
      real(r8) :: pcut    ! ptan threshold for high tangent heights
                                        
      ! retrieval work arrays
      real(r8) :: sensitivity                         ! sensitivity with slope and correction
      real(r8) :: teff                                ! effective optical depth
      real(r8) :: trans                                ! transmission function
      real(r8) :: y      ! measurement array
      real(r8) :: sy     ! variance of y
      integer :: n1
      logical :: doMaf      ! array for MAF flag
      real(r8), dimension(:), allocatable :: tmp1, tmp2      ! working array
      real(rm), dimension(:,:), allocatable :: A      ! working array
      real(r8), dimension(:), allocatable :: dx       ! working array for x
      real(r8), dimension(:), allocatable :: x        ! s grid array
      real(r8), dimension(:), allocatable :: x0       ! A priori of x
      real(r8), dimension(:), allocatable :: sx0       ! variance of A priori
      real(r8), dimension(:), allocatable :: xext       ! extinction along los
      real(r8), dimension(:,:), allocatable :: sx       ! variance of x
      real(r8), dimension(:,:), allocatable :: C    ! for problem y=Kx, c=K^t#Sy^-1#K
                                                      ! last dimension is for mif
                                                      ! first two are (chan, s)

      ! use this for testing
      pcut = -2.5
      
      ! find how many MAFs        
      nMAFs = chunk%lastMAFIndex-chunk%firstMAFIndex + 1
      
      if (size(configIndices) > 1 .or. &
        & size(configDatabase(configIndices(1))%signals) > 1) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Only one signal is allowed in high cloud retrieval' )

      ! get signal information for this model. Note: allow only 1 signal 
        signal = configDatabase(configIndices(1))%signals(1)
      
        call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows',ModuleName )
      ! create FwdModelOut1 vector by cloning measurements
        call cloneVector ( FwdModelOut1, measurements, vectorNameText='_covarianceDiag' )          

      ! find which channel is used
        ich = 0
        nFreqs = size (signal%frequencies)
        do k=1,nFreqs
              if(signal%channels(k)) ich=k
        end do

      !memorize the initial model configuration
!        cloudysky = configDatabase(configIndices(1))%cloud_width
        cloudysky = configDatabase(configIndices(1))%i_saturation

      ! create covarianceDiag array
        call cloneVector ( covarianceDiag, state, vectorNameText='_covarianceDiag' )
      ! get the inverted diagnonal elements of covariance of apriori
        call getDiagonal ( covariance%m, covarianceDiag )
        
      ! find about the dimension of the state vector        
         xExtPtr => GetVectorQuantityByType ( state, quantityType=l_cloudextinction)
         nInst = xExtPtr%template%noInstances
         nz = xExtPtr%template%noSurfs
        
      ! allocate C, y, x matrices
          allocate(A(nz*nInst,nz*nInst),C(nz*nInst,nz*nInst))
          allocate(x(nz*nInst),x0(nz*nInst),sx0(nz*nInst),dx(nz*nInst))
          allocate(tmp1(nz),tmp2(nz))
          A = 0._rm
          C = 0._r8
          y = 0._r8
          sy = 0._r8
          dx = 0._r8
            
          ! get cloud radiance measurements for this signal
          Tcir => GetVectorQuantityByType ( Measurements,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )
            
          ! get cloud radiances error for this signal
          Terr => GetVectorQuantityByType ( MeasurementSD,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )

          ! get pointers of x covariance for the retrieval
          xExtVar => GetVectorQuantityByType (covarianceDiag, &
               & quantityType=l_cloudextinction)
      
          ! get pointers of output SD for the retrieval
          outExtSD => GetVectorQuantityByType (outputSD, &
               & quantityType=l_cloudextinction)
      
          ! get cloud radiance measurements for this signal
          ConstrainTcir => GetVectorQuantityByType ( fwdModelExtra,     &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )
          
          ConstrainTcir%values = Tcir%values
          
      ! Loop over MAFs
      do maf=1, nMAFs
      fmStat%maf = maf
                        
          ! to save time, skip cloud sensitivity calculation if there is no cloud
          ! overwrite model configuration
          doMaf = .false.
          if(sum(ConstrainTcir%values(:,maf)) .ne. 0._r8) doMaf = .true.
          
!          configDatabase(configIndices(1))%cloud_width=0
!          if(doMaf) configDatabase(configIndices(1))%cloud_width=cloudysky
          configDatabase(configIndices(1))%i_saturation=0
          if(doMaf) configDatabase(configIndices(1))%i_saturation=cloudysky

      print*,'begin cloud retrieval maf= ',maf,' chunk size=',nMAFs,'type=',doMaf
          fmStat%rows = .false.
          call forwardModel ( configDatabase(configIndices(1)), &
            & state, fwdModelExtra, FwdModelOut1, fmw, fmStat, jacobian )
            
          ! get clear sky radiances from forward model for this signal
          Tb0 => GetVectorQuantityByType ( fwdModelOut1,                 &
               & quantityType=l_radiance,                                      &
               & signal=signal%index, sideband=signal%sideband )

          ! get ptan.  Warning: all signals must be in the same module
          ptan => GetVectorQuantityByType ( fwdModelExtra,      &
               & quantityType=l_ptan, instrumentModule = &
               & Tb0%template%instrumentModule)

          ! get earthradius from forward model
!          Re => GetVectorQuantityByType ( fwdModelExtra,      &
!               & quantityType=l_earthradius)
          
          ! get sensitivity from forward model for this signal
          Slope => GetVectorQuantityByType ( fwdModelOut1,      &
               & quantityType=l_cloudRADSensitivity,                           &
               & signal=signal%index, sideband=signal%sideband )
          
          ! get rowBlock and colBlock for this model
          rowJBlock = FindBlock (jacobian%row, Tb0%index, maf)
          colJBlock = FindBlock ( Jacobian%col, xExtPtr%index, maf )
          jBlock => jacobian%block(rowJblock,colJblock)
               
          do mif=1,ptan%template%noSurfs
          if(ptan%values(mif,maf) > pcut) then
             
             sensitivity = 0._r8
             if(doMaf) then
                ! find more accurate sensitivity
                do j=1,nz
                  tmp1(j) = 1._r8/nz*(j-1._r8)
                  tmp2(j) = slope%values(ich+nFreqs*(mif-1),maf)* &
                     & (1._r8 - 0.2_r8* tmp1(j)**5) ! correction term (see ATBD)
!		            if(j > 1) then
!			            if(tmp2(j) > tmp2(j-1)) n1=j
!		            end if
                end do
!		          if(n1 < nz/2) print*,'sensitivity too low for mif=',mif
                     
                y = Tcir%values(ich+nFreqs*(mif-1),maf)

                ! interpolate to get initial guess of teff
                teff = 0._r8
                if(y < tmp2(1)) teff=0._r8
                if(y > tmp2(nz)) teff=1._r8
                do j=1,nz-1
                  if(y > tmp2(j) .and. y < tmp2(j+1)) then
                    teff=(tmp1(j)*(y-tmp2(j))+tmp1(j+1)*(tmp2(j+1)-y)) &
                      & /(tmp2(j+1)-tmp2(j))
                  end if
                end do

                ! solve the relation one more time to refine teff and sensitivity
                sensitivity = slope%values(ich+nFreqs*(mif-1),maf)* &
                   & (1._r8 - 0.2_r8* teff**5) ! correction term (see ATBD)
            end if  ! when sensitivity is calculated

                ! in case we have small sensitivity or no cloud
                if(sensitivity < 0._r8) cycle
                if(abs(sensitivity) < 1._r8) sensitivity = 130._r8
                     
                teff = y/sensitivity
                if(teff > 1._r8) teff = 1._r8
                           
                ! convert cloud radiance to effective optical depth
                y = teff
                sy = Terr%values(ich+nFreqs*(mif-1),maf)**2/sensitivity**2
                     
                do i=1,nz
                  do j=1,nInst
                    dx(i+(j-1)*nz) = dx(i+(j-1)*nz) + jBlock%values(mif,i+(j-1)*nz)*y/sy
                  end do
                end do
                do i=1,nz
                 do j=1,nInst
                  do i1=1,nz
                   do j1=1,nInst
                   C(i+(j-1)*nz,i1+(j1-1)*nz) = C(i+(j-1)*nz,i1+(j1-1)*nz) + &
                    & jBlock%values(mif,i+(j-1)*nz)*jBlock%values(mif,i1+(j1-1)*nz)/sy
                   end do
                  end do
                 end do
                end do
          end if
          end do   ! mif

      call clearMatrix ( jacobian )           ! free the space
      end do ! end of mafs
      
    ! give back the model config value
!      configDatabase(configIndices(1))%cloud_width=cloudysky
      configDatabase(configIndices(1))%i_saturation=cloudysky

      ! start inversion
      x0 = 0._r8        ! A priori
      x = 0._r8
      A = C
      do i=1,nz
         do j=1,nInst
         sx0(i+(j-1)*nz) = xExtVar%values(i,j)   ! sx has already been inverted
         A(i+(j-1)*nz,i+(j-1)*nz)=A(i+(j-1)*nz,i+(j-1)*nz) + xExtVar%values(i,j) 
         end do
      end do
         
      Call MatrixInversion(A)
               
      ! output estimated SD 
      do i=1,nz
         do j=1,nInst
            outExtSD%values(i,j) = A(i+(j-1)*nz,i+(j-1)*nz)
         end do
      end do

      ! compute the least-squared solution
      do i=1,nz
         do j=1,nInst
         x(i+(j-1)*nz) = sum(reshape(A(i+(j-1)*nz,:),(/nz*nInst/))*(x0*sx0 + dx))
         end do
      end do
      
      ! output retrieval results
      xExtPtr%values = reshape(x,(/nz,nInst/))
      
      ! deallocate arrays and free memory
      deallocate(A,C,tmp1,tmp2,dx,x,x0,sx0)

    ! ExtSD is defined as MLS contribution only
           badValue = outExtSD%template%badValue
!      outExtSD%values = 1._r8/outExtSD%values - xExtVar%values !xExtVar is an inverted variance
!       where(outExtSD%values > 0._r8) outExtSD%values = 1._r8/sqrt(outExtSD%values)
!       where(outExtSD%values .le. 0._r8) outExtSD%values = badValue

   ! clean up
      call destroyVectorInfo ( FwdModelOut1 )
      call destroyVectorInfo ( CovarianceDiag )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      
    end subroutine HighCloudRetrieval

! ------------------------------------------  LowCloudRetrieval  -----
    subroutine LowCloudRetrieval
      
      ! Local Variables
      type (ForwardModelIntermediate_T) :: Fmw     ! Work space for forward model
      type (VectorValue_T), pointer :: xLosPtr        ! pointer of l_lostransfunc quantity
      type (VectorValue_T), pointer :: xExtPtr        ! pointer of l_cloudExtinction quantity
      type (VectorValue_T), pointer :: xLosVar        ! variance of apriori
      type (VectorValue_T), pointer :: xExtVar        ! variance of apriori
      type (VectorValue_T), pointer :: outLosSD        ! SD of output
      type (VectorValue_T), pointer :: outExtSD        ! SD of output
      type (VectorValue_T), pointer :: ConstrainTcir        ! for model cloud top indicator
      type (VectorValue_T), pointer :: Tcir        ! cloud-induced radiance
      type (VectorValue_T), pointer :: Terr        ! cloud-induced radiance SD
      type (VectorValue_T), pointer :: PTAN        ! Tgt pressure
      type (VectorValue_T), pointer :: Re        ! Earth Radius
      type (VectorValue_T), pointer :: Tb0         ! model clear sky radiance (100%RH)
      type (VectorValue_T), pointer :: Slope       ! sensitivity slope to convert cloud
                                                   ! radiance to optical depth
                                          
      integer :: i,j,k,imodel,mif,maf,isignal          ! Loop subscripts
      integer :: ich, ich0          ! channel used
      integer :: cloudysky(size(configIndices))   ! cloudysky index from Model Configuration

      integer :: nFreqs0      ! number of frequencies in the signal for cloud top estimation
      integer :: nSignal      ! number of signals in a band  
      integer :: nSgrid      ! number of S grids  
      integer :: nChans      ! total number of channels used in retrieval
      integer :: nMifs       ! number of total mifs
      integer :: ndoMifs       ! number of used mifs
      real(r8) :: p_lowcut    ! ptan threshold for low tangent heights
      real(r8) :: scale
                                        
      ! retrieval work arrays
      real(r8) :: sensitivity                         ! sensitivity with slope and correction
      real(r8) :: teff                                ! effective optical depth
      real(r8) :: trans                                ! transmission function
      integer :: itop                                ! cloud top tangt pressure index
      real(r8) :: zt
      logical, dimension(:), allocatable :: doMaf      ! array for MAF flag
      real(rm), dimension(:,:), allocatable :: A      ! working array
      real(r8), dimension(:,:), allocatable :: dx       ! working array for x
      real(r8), dimension(:,:), allocatable :: y      ! measurement array
      real(r8), dimension(:,:), allocatable :: sy     ! variance of y
      real(r8), dimension(:), allocatable :: Slevel        ! s grid
      real(r8), dimension(:), allocatable :: x        ! s grid array
      real(r8), dimension(:), allocatable :: x0       ! A priori of x
      real(r8), dimension(:), allocatable :: sx0       ! variance of A priori
      real(r8), dimension(:), allocatable :: xext       ! extinction along los
      real(r8), dimension(:,:), allocatable :: sx       ! variance of x
      real(r8), dimension(:,:,:), allocatable :: C    ! for problem y=Kx, c=K^t#Sy^-1#K
                                                      ! last dimension is for mif
                                                      ! first two are (chan, s)
! ----- executable ----
      
      p_lowcut = -2.5
      ! find how many MAFs        
      nMAFs = chunk%lastMAFIndex-chunk%firstMAFIndex + 1
            
      allocate(doMaf(nMAFs))
               
      call allocate_test ( fmStat%rows, jacobian%row%nb, 'fmStat%rows', ModuleName )
      ! create FwdModelOut1 vector by cloning measurements
        call cloneVector ( FwdModelOut1, measurements, vectorNameText='_covarianceDiag' )
               
      ! Measurements are cloud radiances and will be converted to
      ! the effective cloud optical depth after the forward model
      ! provides the sensitivities.
        
      ! Get  l_lostransfunc quantity from state vector
      ! State vector should contains only one quantity of l_lostransfunc type 
      ! (a los quantity),which is the increment of cloud transmission function,
      ! and it is going to be retrieved here.
      
        xLosPtr => GetVectorQuantityByType ( state, quantityType=l_lostransfunc)
        xExtPtr => GetVectorQuantityByType ( state, quantityType=l_cloudextinction)
      
      ! Get S grid dimensions
        nSgrid=xLosPtr%template%noChans    ! number of S grids
        allocate(Slevel(nSgrid))
        sLevel = xLosPtr%template%frequencies      ! sLevel has unit of km here                 
        
        ! get Tcir from extraState vector for constraining cloud top
        ConstrainTcir => GetVectorQuantityByType ( fwdModelExtra,     &
         & quantityType=l_cloudInducedRadiance)
        
      ! find how many channels total in all models
        doMaf = .false.
        nChans = 0
        do imodel = 1, size(configIndices)
          !memorize the initial model configuration
           cloudysky(imodel) = configDatabase(configIndices(imodel))%i_saturation
        do isignal = 1, size(configDatabase(configIndices(imodel))%signals)
          nChans = nChans + &
          & count(configDatabase(configIndices(imodel))%signals(isignal)%channels)
              ! get signal information for this model. Note: allow only 1 signal 
              signal = configDatabase(configIndices(imodel))%signals(isignal)
              ! get cloud radiance measurements for this signal
              Tcir => GetVectorQuantityByType ( Measurements,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband)
              ! determine MAF flag (i.e., whether to call Forward Model for the MAF)
              if(any(iand(ichar(Tcir%mask(:, maf)), m_cloud) == 1)) doMaf(maf) = .true.

        end do ! band signals
        end do ! configIndices or models

      ! find how many mifs from the first quantity        
        nMifs = measurements%quantities(1)%template%noSurfs
               
      ! create covarianceDiag array
        call cloneVector ( covarianceDiag, state, vectorNameText='_covarianceDiag' )
      ! get the inverted diagnonal elements of covariance of apriori
        call getDiagonal ( covariance%m, covarianceDiag )
        
      ! get pointers of x covariance for the apriori
        xLosVar => GetVectorQuantityByType(covarianceDiag,quantityType=l_lostransfunc)
        xExtVar => GetVectorQuantityByType(covarianceDiag,quantityType=l_cloudextinction)
      
      ! get pointers of output SD for the retrieval
        outLosSD => GetVectorQuantityByType(outputSD,quantityType=l_lostransfunc)
        outExtSD => GetVectorQuantityByType(outputSD,quantityType=l_cloudextinction)
        outLosSD%values = 0._r8
        outExtSD%values = 0._r8
      
      ! allocate C, y, x matrices
          allocate(A(nSgrid,nSgrid),C(nSgrid,nSgrid,nMifs))
          allocate(y(nChans,nMifs),sy(nChans,nMifs))
          allocate(x(nSgrid),x0(nSgrid),sx0(nSgrid),xext(nSgrid))
          allocate(dx(nSgrid,nMifs),sx(nSgrid,nMifs))

      ! Loop over MAFs
        do maf =1,nMAFs
        fmStat%maf = maf
        print*,'begin cloud retrieval maf= ',maf,' chunk size=',nMAFs, 'type= ',&
          & configDatabase(configIndices(1))%i_saturation
                        
          A = 0._rm
          C = 0._r8
          y = 0._r8
          sy = 0._r8
          dx = 0._r8
            
          ich = 0
          do imodel = 1, size(configIndices)
            fmStat%rows = .false.
            nSignal = size(configDatabase(configIndices(imodel))%signals)
            ! to save time, skip cloud sensitivity calculation if there is no cloud
            ! overwrite model configuration
            configDatabase(configIndices(imodel))%i_saturation=l_clear
            if(doMaf(maf)) configDatabase(configIndices(imodel))%i_saturation=cloudysky(imodel)
                       
            ! use the cloud radiance in last signal for cloud top indicator 
            signal = configDatabase(configIndices(imodel))%signals(nSignal)

              Tcir => GetVectorQuantityByType ( Measurements,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )
            
!              ConstrainTcir%values = Tcir%values
 
              nFreqs0 = size(signal%frequencies)
               do k=1,nFreqs0
                  if(signal%channels(k)) ich0 = k
               end do
               
            ! estimate cloud top from the Tcir in previous scans
            ! (in simulation, now use the current scan)
              itop = 0
              do mif = 1, nMifs
               if(ConstrainTcir%values(ich0+(mif-1)*nFreqs0,maf) .ne. 0._r8) &
                & itop = mif
              end do
              
            ! call full cloud model for Jacobian and Sensitivity  
             call forwardModel ( configDatabase(configIndices(imodel)), &
                & state, fwdModelExtra, FwdModelOut1, fmw, fmStat, jacobian )
            
            do isignal = 1, nSignal
              ! get signal information for this model. Note: allow only 1 signal 
              signal = configDatabase(configIndices(imodel))%signals(isignal)

              ! get clear sky radiances from forward model for this signal
              Tb0 => GetVectorQuantityByType ( fwdModelOut1,                 &
               & quantityType=l_radiance,                                      &
               & signal=signal%index, sideband=signal%sideband )

              ! get ptan.  Warning: all signals must be in the same module
              ptan => GetVectorQuantityByType ( fwdModelExtra,      &
               & quantityType=l_ptan, instrumentModule = &
               & Tb0%template%instrumentModule)

              ! get earthradius from forward model
              Re => GetVectorQuantityByType ( fwdModelExtra,      &
               & quantityType=l_earthradius)
          
              ! get sensitivity from forward model for this signal
              Slope => GetVectorQuantityByType ( fwdModelOut1,      &
               & quantityType=l_cloudRADSensitivity,                           &
               & signal=signal%index, sideband=signal%sideband )
          
              ! get cloud radiance measurements for this signal
              Tcir => GetVectorQuantityByType ( Measurements,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )
            
              ! get cloud radiances error for this signal
              Terr => GetVectorQuantityByType ( MeasurementSD,       &
               & quantityType=l_cloudInducedRadiance,             &
               & signal=signal%index, sideband=signal%sideband )

              ! get rowBlock and colBlock for this model
              rowJBlock = FindBlock (jacobian%row, Tb0%index, maf)
              colJBlock = FindBlock ( Jacobian%col, xLosPtr%index, maf )
                
              jBlock => jacobian%block(rowJblock,colJblock)
               
              nFreqs = size (signal%frequencies)
               do k=1,nFreqs
               if(signal%channels(k)) then
               ich = ich + 1
                  do mif=1,nMifs
                  if(ptan%values(mif,maf) < p_lowcut) then
                    
                    sensitivity = 0._r8
                    if(doMaf(maf)) then
                     ! find more accurate sensitivity. we first borrow 
                     ! x, x0 arrays to establish the Tcir-teff relation:
                     ! x-> Tcir; x0->teff;
                     do j=1,nSgrid
                     x0(j) = 1._r8/nSgrid*(j-1._r8)
                     x(j) = slope%values(k+nFreqs*(mif-1),maf)* &
                        & (1._r8 + 0.46_r8* x0(j)**6) ! correction term (see ATBD)
                     end do
                     
                     y(ich,mif) = Tcir%values(k+nFreqs*(mif-1),maf)

                     ! interpolate to get initial guess of teff
                     teff = 0._r8
                     if(y(ich,mif) < x(1)) teff=0._r8
                     if(y(ich,mif) > x(nSgrid)) teff=1._r8
                     do j=1,nSgrid-1
                        if(y(ich,mif) > x(j) .and. y(ich,mif) < x(j+1)) then
                        teff=(x0(j)*(y(ich,mif)-x(j))+x0(j+1)*(x(j+1)-y(ich,mif))) &
                          & /(x(j+1)-x(j))
                        end if
                     end do

                     ! solve the relation one more time to refine teff and sensitivity
                     sensitivity = slope%values(k+nFreqs*(mif-1),maf)* &
                           & (1._r8 + 0.46_r8* teff**6) ! correction term (see ATBD)
                    end if ! when sensitivity is calculated

                     ! in case we have small sensitivity or no cloud
                     if(abs(sensitivity) < 1._r8) sensitivity = -105._r8
                     
                     teff = y(ich,mif)/sensitivity
                     if(teff > 1._r8) teff = 1._r8
                           
                     ! convert cloud radiance to effective optical depth
                     
                     y(ich,mif) = teff
                     sy(ich,mif) = (Terr%values(k+nFreqs*(mif-1),maf)/sensitivity)**2
                     
                     do i=1,nSgrid
                        dx(i,mif) = dx(i,mif) + jBlock%values(k,i+(mif-1)*nSgrid)* &
                           & y(ich,mif)/sy(ich,mif)
                     do j=1,nSgrid
                        C(i,j,mif) = C(i,j,mif) + jBlock%values(k,i+(mif-1)*nSgrid)* &
                           & jBlock%values(k,j+(mif-1)*nSgrid)/sy(ich,mif)
                           
                     end do
                     end do
                  end if
                  end do   ! mif

               end if
               end do      ! end of frequency k

            end do         ! end of band signals
         end do         ! end of imodel

         ! check if Jacobian rows are consistent with Signal rows
         if(ich /= nChans) call MLSMessage ( MLSMSG_Warning, ModuleName, &
           & 'inconsistent channels between Jacobian and Signal' )

           sx = 1.e8_r8         ! sx is the inversd variance of a priori
           x = 0._r8
           x0 = xLosVar%template%frequencies
           do mif=1,nMifs
             do i=1,nSgrid
                !xLosVar is the inverted variance
                sx(i,mif) = xLosVar%values(i+nSgrid*(mif-1),maf)  ! xLosVar has been inverted
             end do
             if(itop .ne. 0) then
             !constrained by the estimated cloud top
                x = x0**2/2._r8/(re%values(1,maf)*0.001_r8 + &
                  & (ptan%values(mif,maf)+3._r8)*16._r8)
                x = x/16._r8 + ptan%values(mif,maf)
                do i = 1,nSgrid
                 if(x(i) > ptan%values(itop,maf)) sx(i,mif) = sx(i,mif)*1.e4
                end do
             end if
           end do
                              
         ! start inversion
           do mif=1,nMifs
            if (ptan%values(mif,maf) < p_lowcut) then
               ! for cloud retrieval: x_star=0, y_star=0, x0=apriori=0

               x0 = 0._r8        ! A priori
               x = 0._r8
               sx0=sx(:,mif)
               
               do j=1,maxJacobians     ! maxJacobians is number of iterations (default 5)
               
                  x0 = x        ! use  the last retrieval as A priori
                                 ! and doubt constraints to the A priori
                  where(x0 .le. 0._r8) x0 = 0._r8

                  A = reshape(C(:,:,mif),(/nSgrid,nSgrid/))
                  do i=1,nSgrid
                   A(i,i)=A(i,i) + sx0(i)       ! sx has been inverted
                  end do
                  Call MatrixInversion(A)
               
                  ! output estimated SD after the first iteration
                  if(associated(xLosPtr) .and. j == 1) then
                   do i=1,nSgrid
                    outLosSD%values(i+(mif-1)*nSgrid,maf) = sqrt(A(i,i))
                   end do
                  end if
                  
                  ! compute the LOS retrieval in the least-squard solution
                  do i=1,nSgrid
                     x(i) = sum( reshape(A(i,:),(/nSgrid/))* &
                     & (x0*sx0 + reshape(dx(:,mif),(/nSgrid/)) ))
                  end do
                  ! reduce covariance on x0 for next iteration, using factor lambda (default 10)
                  sx0 = sx0*initLambda          ! sx0 is the inversed variance

               end do  ! end of iteration
               
               ! Now, x is the los Transmission increment
               ! compute cloud extinction from los Transmission increments
               xext = 0._r8      ! temporary storage for LOS extinction
               trans = 0._r8
               do i=nSgrid,2,-1
                 trans = trans + x(i)
                 ! protect from blowing up
                 if(1._r8-trans > 0.004) then
                  scale = (1._r8-trans)*(slevel(i)-slevel(i-1))      ! see ATBD
                                 xext(i)=x(i)/scale
                  ! and standard deviation
                  outLosSD%values(i+(mif-1)*nSgrid,maf) = &     !neglect higher orders
                     & outLosSD%values(i+(mif-1)*nSgrid,maf)/scale
                  xLosVar%values(i+(mif-1)*nSgrid,maf) = &     !xLosVar is the inversed var
                     & xLosVar%values(i+(mif-1)*nSgrid,maf)*scale*scale
                 end if
               end do
               
               ! output los extinction to state vector
               if(associated(xLosPtr)) then
                do i=1,nSgrid
                 xLosPtr%values(i+(mif-1)*nSgrid,maf) = xext(i)
                end do
               end if               
               
            end if
           end do  ! mif
      call clearMatrix ( jacobian )           ! free the space
      
      end do ! end of mafs
      
    ! give back the model config value
      do imodel = 1,size(configIndices)
        configDatabase(configIndices(imodel))%i_saturation=cloudysky(imodel)
      end do
      
    ! deallocate arrays and free memory
      deallocate(A,C,y,sy,dx,x,x0,sx0,xext,sx)

      xExtPtr%values = 0._r8
      outLosSD%values = outLosSD%values**2
      xLosVar%values = 1._r8/xLosVar%values  ! xLosVar is the inverted variance
      
   	call LOS2Grid(xExtPtr,outExtSD,xExtvar,xLosPtr,outLosSD,xLosVar,Ptan,Re,p_lowcut)

    ! output SD
      outLosSD%values = sqrt(outLosSD%values)
    ! ExtSD is defined as MLS contribution only
      badValue = outExtSD%template%badValue
      outExtSD%values = 1._r8/outExtSD%values - 1._r8/xExtVar%values
       where(outExtSD%values > 0._r8) outExtSD%values = 1._r8/sqrt(outExtSD%values)
       where(outExtSD%values .le. 0._r8) outExtSD%values = badValue
    ! Output SD
      xExtVar%values = sqrt(abs(xExtVar%values))
      
    ! overwrite input covariance
      call clearMatrix(covariance%m)     ! this will initialize it to M_Absent
      call updateDiagonal ( covariance, CovarianceDiag)

    ! clean up
      call destroyVectorInfo ( FwdModelOut1 )
      call destroyVectorInfo ( CovarianceDiag )
      call deallocate_test ( fmStat%rows, 'FmStat%rows', moduleName )
      deallocate(Slevel, doMaf)
      
    end subroutine LowCloudRetrieval

  !=============================== LOS2Grid ===================================
  subroutine LOS2Grid( Qty, vQty, vQtyApr, Los, vLos, vLosApr, Ptan, Re, Pcut)

    ! This is to fill a l2gp type of quantity with a los grid type of quantity.
    ! The los quantity is a vector quantity that has dimension of (s, mif, maf),
    ! where s is the path along los.
    !
    ! Linear interpolation is used to fill l2gp grids and unfilled grids are
    ! marked with the baddata flag (-999.)

    use UNITS
    use MLSNumerics, only: InterpolateValues
    use VectorsModule, only: VectorValue_T

    ! Dummy arguments
    type (VectorValue_T), intent(in) :: LOS ! LOS quantity
    type (VectorValue_T), intent(in) :: vLOS ! variance of LOS
    type (VectorValue_T), intent(in) :: vLosApr ! variance of LOS Apriori
    type (VectorValue_T), intent(in) :: Ptan ! tangent pressure
    type (VectorValue_T), intent(in) :: Re ! Earth's radius
    type (VectorValue_T), INTENT(INOUT) :: Qty ! Quantity to fill
    type (VectorValue_T), INTENT(INOUT) :: vQty ! variance of Quantity
    type (VectorValue_T), INTENT(INOUT) :: vQtyApr ! variance of Quantity Aprori

    ! Local variables
    integer :: i, j, maf, mif                ! Loop counter
    integer :: maxZ, minZ                    ! pressure range indices of sGrid
    integer :: noMAFs,noMIFs,noDepths
    real (r8) :: pcut
    real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: cnta
    real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: cnt
    real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: out
    real (r8), dimension(qty%template%noSurfs) :: outZeta, phi_out, beta_out, weight, weighta
    real (r8), dimension(los%template%noChans) :: x_in, y_in, sLevel
    real (r8), dimension(los%template%noSurfs) :: zt

    if ( qty%template%verticalCoordinate == l_pressure ) then
        outZeta = -log10 ( qty%template%surfs(:,1) )
    else
        outZeta = qty%template%surfs(:,1)
    end if

    noMAFs=los%template%noInstances
    noMIFs=los%template%noSurfs
! Now, we use frequency coordinate as sGrid along the path
    noDepths=los%template%noChans
    sLevel = los%template%frequencies
   
! initialize quantity
   Qty%values=Qty%template%badValue
   vQty%values = vQty%template%badValue
   vQtyApr%values = vQtyApr%template%badValue
   cnta=0._r8
   cnt=0._r8
   out=0._r8
   
    do maf=1,noMAFs  
      zt = ptan%values(:,maf)   ! noChans=1 for ptan
      zt = (zt+3.)*16.                      ! converted to height in km
      do mif=1,noMIFs
      if (ptan%values(mif,maf) .gt. pcut) cycle
      ! find altitude of each s grid
      x_in = sLevel**2/2./(re%values(1,maf)*0.001_r8 + zt(mif))
      ! converted to zeta
      x_in = x_in/16. + ptan%values(mif,maf)
      ! find minimum and maximum pressures indices in sGrid
        do i = 2,qty%template%noSurfs-1
        if (ptan%values(mif,maf) < (outZeta(i)+outZeta(i+1))/2. .and. &
          & ptan%values(mif,maf) > (outZeta(i)+outZeta(i-1))/2.) &
          & minZ = i
        end do
        if (ptan%values(mif,maf) < (outZeta(1)+outZeta(2))/2.) minZ=1
        if (ptan%values(mif,maf) > outZeta(qty%template%noSurfs)) cycle ! goto next mif
        
        do i = 2,qty%template%noSurfs-1
        if (x_in(noDepths) < (outZeta(i)+outZeta(i+1))/2. .and. &
          & x_in(noDepths) > (outZeta(i)+outZeta(i-1))/2.) &
          & maxZ = i
        end do
        if (x_in(noDepths) < (outZeta(1)+outZeta(2))/2.) cycle    ! goto next mif
        if (x_in(noDepths) > outZeta(qty%template%noSurfs)) maxZ=qty%template%noSurfs

      ! get phi along path for each mif (phi is in degree)
      y_in = los%template%phi(mif,maf) &
        & - atan(sLevel/(re%values(1,maf)*0.001_r8 + zt(mif)))*180._r8/Pi
       ! interpolate phi onto standard vertical grids     
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),phi_out(minZ:maxZ), &
           & method='Linear')
       ! interpolate quantity to standard vertical grids      
        y_in = Los%values((1+(mif-1)*noDepths):mif*noDepths,maf)
        beta_out = 0._r8
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),beta_out(minZ:maxZ), &
           & method='Linear')
        y_in = vLos%values((1+(mif-1)*noDepths):mif*noDepths,maf)
        weight = 0._r8
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),weight(minZ:maxZ), &
           & method='Linear')
        weight = weight*(maxZ-minZ+1._r8)/size(y_in)  !inflated after interpolation
        y_in = vLosApr%values((1+(mif-1)*noDepths):mif*noDepths,maf)
        weighta = 0._r8
        call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),weighta(minZ:maxZ), &
           & method='Linear')
        weighta = weighta*(maxZ-minZ+1._r8)/size(y_in)  !inflated after interpolation
        
       ! interpolate quantity to standard phi grids
        do i=minZ,maxZ  
          do j = 2, qty%template%noInstances-1
          if(phi_out(i) .lt. &     
            & (qty%template%phi(1,j)+qty%template%phi(1,j+1))/2._r8 &
            & .and. phi_out(i) .ge. &  
            & (qty%template%phi(1,j-1)+qty%template%phi(1,j))/2._r8 ) then
            out(i,j)=out(i,j) + beta_out(i)/weight(i)    ! weighted by variance
            cnt(i,j)=cnt(i,j)+1._r8/weight(i)       !  counter
            cnta(i,j)=cnta(i,j)+1._r8/weighta(i)       !  counter
          end if
          end do
        end do
      end do                            ! End surface loop
    end do                              ! End instance loop
    ! average all non-zero bins
    where (cnt > 0._r8) Qty%values = out/cnt
    where (cnt > 0._r8) cnt = 1._r8/cnt
    where (cnta > 0._r8) cnta = 1._r8/cnta
    
    where (cnta > 0._r8) vQtyApr%values = cnta
    where (cnt > 0._r8) vQty%values = cnt
!    where (cnt > 0._r8 .and. cnt > 0.5_r8*cnta) vQty%values = -cnt  ! assign negative

  end subroutine LOS2Grid    

 end subroutine CloudRetrieval
end module CloudRetrievalModule
! $Log$
! Revision 2.4  2003/05/14 03:55:25  dwu
! tidy up
!
! Revision 2.3  2003/05/13 22:25:52  dwu
! changes in lowcloudretrieval
!
! Revision 2.2  2003/05/13 20:48:38  dwu
! fix a bug
!
! Revision 2.1  2003/05/13 20:38:55  dwu
! change to pass chunk in
!
