module NO_CONV_AT_ALL_M
  use MLSCommon, only: I4, R4, R8
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  use D_LINTRP_M, only: LINTRP
  use D_CSPLINE_M, only: CSPLINE
  use DCSPLINE_DER_M, only: CSPLINE_DER
  use dump_0,only:dump
  use VectorsModule, only: Vector_T, VectorValue_T, GETVECTORQUANTITYBYTYPE
  use ForwardModelConfig, only: ForwardModelConfig_T
  use Intrinsic, only: L_VMR
  use String_Table, only: GET_STRING
  use MatrixModule_0, only: M_BANDED, M_FULL, DUMP
  use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T, DUMP
  implicit NONE
  private
  public :: NO_CONV_AT_ALL
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
! This subroutine transfers the derivatives over from the internal
! convolution grid to the users specified points. This module uses
! cubic spline interpolation to do the job.
!
  Subroutine no_conv_at_all (forwardModelConfig, forwardModelIn, maf, channel, &
    & windowStart, windowFinish, temp, ptan, radiance, &
    & tan_press,i_raw,k_temp,k_atmos, i_star_all, Jacobian )
!
    type (ForwardModelConfig_T) :: FORWARDMODELCONFIG
    type (Vector_T), intent(in) :: FORWARDMODELIN
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    integer, intent(in) :: WINDOWSTART
    integer, intent(in) :: WINDOWFINISH
    type (VectorValue_T), intent(in) :: TEMP
    type (VectorValue_T), intent(in) :: PTAN
    type (VectorValue_T), intent(in) :: RADIANCE

    real(r8), intent(IN) :: I_RAW(:),TAN_PRESS(:)

    Real(r4) :: k_temp(:,:,:)                    ! (Nptg,mxco,mnp)
    Real(r4) :: k_atmos(:,:,:,:)                 ! (Nptg,mxco,mnp,Nsps)

    type (Matrix_T), intent(inout), optional :: Jacobian
!
! -----     Output Variables   ----------------------------------------
!
    real(r8), intent(OUT) :: I_STAR_ALL(:)
!
! -----     Local Variables     ----------------------------------------
!
    integer:: no_t,no_tan_hts,no_phi_t
    type (VectorValue_T), pointer :: F  ! vmr quantity
    integer :: nz
    integer :: N, I, IS, J, K, NF, SV_I, Spectag
    integer :: phiWindow
    integer :: row, col                 ! Indices
    integer :: ptg                      ! Index
!
    real(r8) :: RAD(ptan%template%noSurfs), SRad(ptan%template%noSurfs)
    real(r8) :: der_all(ptan%template%noSurfs)
!
    Character(LEN=01) :: CA
!
! -----  Begin the code  -----------------------------------------
!
    no_t = temp%template%noSurfs
    no_phi_t = temp%template%noInstances
    no_tan_hts = size(tan_press)
!
! Compute the ratio of the strengths
!
! This subroutine is called by channel
!
    k = no_tan_hts
    j = ptan%template%noSurfs
    if (present(Jacobian)) then
      Call Cspline_der(tan_press,Ptan%values(:,maf),i_raw,i_star_all,der_all,k,j)
      row = FindBlock ( Jacobian%row, radiance%index, maf )
      col = FindBlock ( Jacobian%col, ptan%index, maf )
      call CreateBlock ( Jacobian, row, col, m_banded, &
        & radiance%template%noSurfs*radiance%template%noChans )
      jacobian%block(row,col)%values = 0.0_r8
      do ptg = 1, j
        jacobian%block(row,col)%values( channel + &
          & radiance%template%noChans*(ptg-1),1 ) = der_all(ptg)
        jacobian%block(row,col)%r1(ptg) = &
          & 1 + radiance%template%noChans*(ptg-1)
        jacobian%block(row,col)%r2(ptg) = &
          & radiance%template%noChans*ptg
      end do
    else
      Call Cspline(tan_press,Ptan%values(:,maf),i_raw,i_star_all,k,j)
    endif
      
!
    if(.not. ANY((/forwardModelConfig%temp_der,forwardModelConfig%atmos_der,&
      & forwardModelConfig%spect_der/))) Return
    row = FindBlock ( Jacobian%row, radiance%index, maf ) ! Tidy up conditions later !???

!
! Now transfer the other fwd_mdl derivatives to the output pointing
! values
!
! ********************* Temperature derivatives ******************
!
! check to determine if derivative is desired for this parameter
!
    if (forwardModelConfig%temp_der) then
!
! Derivatives needed continue to process
!
      Rad(1:) = 0.0
      phiWindow = windowFinish-windowStart+1
      do nf = 1, phiWindow
        col = FindBlock ( Jacobian%col, temp%index, nf+windowStart-1 )
        call CreateBlock ( Jacobian, row, col, m_full )

        do sv_i = 1, no_t
          Rad(1:k) = k_temp(1:k,sv_i,nf)
          Call Cspline(tan_press,Ptan%values(:,maf),Rad,SRad,k,j)
          do ptg = 1,j
            jacobian%block(row,col)%values( channel+&
              & radiance%template%noChans*(ptg-1),sv_i) = Srad(ptg)
          end do
        end do

      end do


    endif

    if(forwardModelConfig%atmos_der) then

      ! ****************** atmospheric derivatives ******************

      do is = 1, size(ForwardModelConfig%molecules)

         f => GetVectorQuantityByType ( forwardModelIn, quantityType=l_vmr, &
          & molecule=forwardModelConfig%molecules(is), noError=.true. )

        if (associated(f)) then

          nz = f%template%noSurfs
          Rad(1:) = 0.0

          ! Derivatives needed continue to process

          do nf = 1, f%template%noInstances
            col = FindBlock ( Jacobian%col, f%index, nf+windowStart-1 )
            call CreateBlock ( Jacobian, row, col, m_full )

            do sv_i = 1, f%template%noSurfs
              Rad(1:k) = k_atmos(1:k,sv_i,nf,is)
              Call Lintrp(tan_press,Ptan%values(:,maf),Rad,SRad,k,j)
              do ptg = 1,j
                jacobian%block(row,col)%values( channel+&
                  & radiance%template%noChans*(ptg-1),sv_i) = Srad(ptg)
              end do
            end do

          end do
!
        endif
!
      end do
!
    endif
!
!     if(spect_der) then
!
! ! ****************** Spectroscopic derivatives ******************
! !
!       do is = 1, n_sps
! !
!         i = spect_atmos(is)
!         if(i < 1) CYCLE
!         if(.not.spectroscopic(i)%DER_CALC(band)) CYCLE
! !
! ! Derivatives needed continue to process
! !
!         Spectag = atmospheric(is)%spectag
! !
!         DO
! !
!           if(spectroscopic(i)%Spectag /= Spectag) EXIT
!           n = spectroscopic(i)%no_phi_values
!           nz = spectroscopic(i)%no_zeta_values
!           CA = spectroscopic(i)%type
!           ki = ki + 1
!           kc = kc + 1
!           k_star_info(kc)%name = spectroscopic(i)%NAME
!           k_star_info(kc)%first_dim_index = ki
!           k_star_info(kc)%no_phi_basis = n
!           k_star_info(kc)%no_zeta_basis = nz
!           k_star_info(kc)%zeta_basis(1:nz) = &
!                   &  spectroscopic(i)%zeta_basis(1:nz)
! !
!           Rad(1:) = 0.0
!           do nf = 1, n
! !
!             do sv_i = 1, nz
! !
!               select case ( CA )
!                 case ( 'W' )
!                   Rad(1:k) = k_spect_dw(1:k,sv_i,nf,i)
!                 case ( 'N' )
!                   Rad(1:k) = k_spect_dn(1:k,sv_i,nf,i)
!                 case ( 'V' )
!                   Rad(1:k) = k_spect_dnu(1:k,sv_i,nf,i)
!               end select
! !
!               Call Lintrp(tan_press,Ptan,Rad,SRad,k,j)
!               k_star_all(ki,sv_i,nf,1:j) = SRad(1:j)
! !
!             end do        ! sv_i loop
! !
!           end do          ! nf loop
! !
!           i = i + 1
!           if(i > 3 * n_sps) EXIT
! !
!         END DO
! !
!       end do
! !
!     endif
!

    Return
!
  End Subroutine NO_CONV_AT_ALL
!
end module NO_CONV_AT_ALL_M
! $Log$
! Revision 1.10  2001/04/19 23:56:52  livesey
! New parameters
!
! Revision 1.9  2001/04/10 10:14:16  zvi
! Fixing bug in convolve routines
!
! Revision 1.8  2001/04/10 02:25:14  livesey
! Tidied up some code
!
! Revision 1.7  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.6  2001/03/29 02:54:49  livesey
! Removed print statements
!
! Revision 1.5  2001/03/29 02:54:29  livesey
! Changed assumed size to assumed shape
!
! Revision 1.4  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.3  2001/03/21 01:10:38  livesey
! Now gets Ptan from vector
!
! Revision 1.2  2001/03/07 23:45:15  zvi
! Adding logical flags fro Temp, Atmos & Spect. derivatives
!
! Revision 1.1  2000/06/21 21:56:14  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
