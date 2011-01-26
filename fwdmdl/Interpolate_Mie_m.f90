! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Interpolate_Mie_m

  implicit NONE
  private
  public :: Interpolate_Mie

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------  Interpolate_Mie  -----

  subroutine Interpolate_Mie ( Frq_Ind, Eta_T_Path, Eta_IWC_a_Path,      &
                             & Eta_IWC_s_Path, Atmos_Der, Temp_Der,      &
                             & Beta_c_e_path, Beta_c_s_path,             &
                             & dBeta_c_e_dIWC_path, dBeta_c_s_dIWC_path, &
                             & dBeta_c_e_dT_path, dBeta_c_s_dT_path )

    ! Interpolate the Mie tables for Frq_Ind to path IWC and Temperature
    use Get_Eta_Matrix_m, only: Eta_D_T, Interpolate_Stru
    use MLSKinds, only: RP
    use Read_Mie_m, only: dBeta_dIWC_c_a, dBeta_dIWC_c_s, &
      & dBeta_dT_c_a, dBeta_dT_c_s, Log_Beta_c_a, Log_Beta_c_s, Log_Mie

    integer, intent(in) :: Frq_Ind     ! Frequency index for Mie tables
    type(eta_d_t), intent(in) :: Eta_T_Path(:), & ! T coeffs
      & Eta_IWC_a_Path(:), Eta_IWC_s_Path(:) ! Absorption, scattering Coeffs
    logical, intent(in) :: Atmos_Der, Temp_Der ! Compute derivatives?
    real(rp), intent(out) :: Beta_c_s_path(:), Beta_c_e_path(:)
    real(rp), intent(out) :: dBeta_c_s_dIWC_path(:), dBeta_c_e_dIWC_path(:)
    real(rp), intent(out) :: dBeta_c_s_dT_path(:), dBeta_c_e_dT_path(:)

    real(rp) :: Beta_c_a_path(size(Beta_c_s_path))
    real(rp) :: dBeta_c_a_dIWC_path(size(Beta_c_s_path))
    real(rp) :: dBeta_c_a_dT_path(size(Beta_c_s_path))

    call log_mie ! Get log_beta_c_a and log_beta_c_s if not done yet

    ! Interpolate Mie beta_c_a and beta_c_s to T and IWC on path
    ! using interpolating coefficients from Mie tables to T and IWC.
    call interpolate_stru ( log_beta_c_a(:,:,frq_ind), &
      & eta_t_path, eta_iwc_a_path, &
      & beta_c_a_path ) ! Actually getting log(beta_c_a_path)
    beta_c_a_path = exp(beta_c_a_path)
    call interpolate_stru ( log_beta_c_s(:,:,frq_ind), &
      & eta_t_path, eta_iwc_s_path, &
      & beta_c_s_path ) ! Actually getting log(beta_c_s_path)
    beta_c_s_path = exp(beta_c_s_path)

    beta_c_e_path = beta_c_a_path + beta_c_s_path

    if ( atmos_der ) then
      call interpolate_stru ( dBeta_dIWC_c_a(:,:,frq_ind), &
        & eta_t_path, eta_iwc_a_path, &
        & dBeta_c_a_dIWC_path )
      call interpolate_stru ( dBeta_dIWC_c_s(:,:,frq_ind), &
        & eta_t_path, eta_iwc_s_path, &
        & dBeta_c_s_dIWC_path )
      dBeta_c_e_dIWC_path = dBeta_c_a_dIWC_path + dBeta_c_s_dIWC_path
    end if

    if ( temp_der ) then
      call interpolate_stru ( dBeta_dT_c_a(:,:,frq_ind), &
        & eta_t_path, eta_iwc_a_path, &
        & dBeta_c_a_dT_path )
      call interpolate_stru ( dBeta_dT_c_s(:,:,frq_ind), &
        & eta_t_path, eta_iwc_s_path, &
        & dBeta_c_s_dT_path )
      dBeta_c_e_dT_path = dBeta_c_a_dT_path + dBeta_c_s_dT_path
    end if

  end subroutine Interpolate_Mie

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Interpolate_Mie_m

! $Log$
! Revision 2.3  2011/01/26 03:04:13  vsnyder
! Different etas for IWC_a and IWC_s, get beta_c_e stuff from beta_c_a+beta_c_s
!
! Revision 2.2  2010/11/05 20:28:34  vsnyder
! Delete unused declarations
!
! Revision 2.1  2010/08/19 19:17:27  vsnyder
! Initial commit
!
