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

  use Sparse_m, only: Sparse_t
  implicit none
  private
  public :: Interpolate_Mie

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------  Interpolate_Mie  -----

  subroutine Interpolate_Mie ( Frq_Ind, Eta_T_IWC_Path, Atmos_Der, Temp_Der, &
                             & Beta_c_e_path, Beta_c_s_path,                 &
                             & dBeta_c_e_dIWC_path, dBeta_c_s_dIWC_path,     &
                             & dBeta_c_e_dT_path, dBeta_c_s_dT_path )

    ! Interpolate the Mie tables for Frq_Ind to path Temperature and IWC
    use MLSKinds, only: RP
    use Read_Mie_m, only: dBeta_dIWC_c_a, dBeta_dIWC_c_s, &
      & dBeta_dT_c_a, dBeta_dT_c_s, Log_Beta_c_a, Log_Beta_c_s, Log_Mie
    use Sparse_Eta_m, only: Sparse_Eta_t

    integer, intent(in) :: Frq_Ind       ! Frequency index for Mie tables
    type(sparse_eta_t) :: Eta_T_IWC_Path ! T x IWC interpolating coeffs
    logical, intent(in) :: Atmos_Der, Temp_Der ! Compute derivatives?
    real(rp), intent(out) :: Beta_c_s_path(:), Beta_c_e_path(:)
    real(rp), intent(out) :: dBeta_c_s_dIWC_path(:), dBeta_c_e_dIWC_path(:)
    real(rp), intent(out) :: dBeta_c_s_dT_path(:), dBeta_c_e_dT_path(:)

    real(rp) :: Beta_c_a_path(size(Beta_c_s_path))
    real(rp) :: dBeta_c_a_dIWC_path(size(Beta_c_s_path))
    real(rp) :: dBeta_c_a_dT_path(size(Beta_c_s_path))
    integer :: NT ! Mie table T x IWC sizes (they're all the same size)
    real(rp), pointer :: Temp_1D(:)

    call log_mie ! Get log_beta_c_a and log_beta_c_s if not done yet

    nt = size(log_beta_c_a(:,:,frq_ind))
    ! Interpolate Mie log_beta_c_a and log_beta_c_s to T and IWC on path
    ! using interpolating coefficients from Mie tables to T and IWC, then
    ! exponentiate to get beta_c_a_path and beta_c_s_path.
    temp_1d(1:nt) => log_beta_c_a(:,:,frq_ind)
    call eta_T_IWC_path%sparse_dot_vec ( temp_1d, beta_c_a_path )
    beta_c_a_path = exp(beta_c_a_path)
    temp_1d(1:nt) => log_beta_c_s(:,:,frq_ind)
    call eta_T_IWC_path%sparse_dot_vec ( temp_1d, beta_c_s_path )
    beta_c_s_path = exp(beta_c_s_path)

    beta_c_e_path = beta_c_a_path + beta_c_s_path

    if ( atmos_der ) then
      temp_1d(1:nt) => dBeta_dIWC_c_a(:,:,frq_ind)
      call eta_T_IWC_path%sparse_dot_vec ( temp_1d, dBeta_c_a_dIWC_path )
      temp_1d(1:nt) => dBeta_dIWC_c_s(:,:,frq_ind)
      call eta_T_IWC_path%sparse_dot_vec ( temp_1d, dBeta_c_s_dIWC_path )
      dBeta_c_e_dIWC_path = dBeta_c_a_dIWC_path + dBeta_c_s_dIWC_path
    end if

    if ( temp_der ) then
      temp_1d(1:nt) => dBeta_dT_c_a(:,:,frq_ind)
      call eta_T_IWC_path%sparse_dot_vec ( temp_1d, dBeta_c_a_dT_path )
      temp_1d(1:nt) => dBeta_dT_c_s(:,:,frq_ind)
      call eta_T_IWC_path%sparse_dot_vec ( temp_1d, dBeta_c_s_dT_path )
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
! Revision 2.6  2018/09/07 17:15:00  pwagner
! Code around a NAG-6.1 compiler bug
!
! Revision 2.5  2018/05/17 02:15:45  vsnyder
! Use sparse instead of dense interpolation
!
! Revision 2.4  2011/07/29 01:57:04  vsnyder
! Only IWC instead of IWC_A and IWC_S
!
! Revision 2.3  2011/01/26 03:04:13  vsnyder
! Different etas for IWC_a and IWC_s, get beta_c_e stuff from beta_c_a+beta_c_s
!
! Revision 2.2  2010/11/05 20:28:34  vsnyder
! Delete unused declarations
!
! Revision 2.1  2010/08/19 19:17:27  vsnyder
! Initial commit
!
