module ZATMOS_DERIV_M
  use MLSCommon, only: I4, R8
  use ForwardModelConfig, only: ForwardModelConfig_T
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use GET_DRAD_NOTDER_M, only: GET_DRAD_NOTDER
  use PATH_ENTITIES_M, only: PATH_DERIVATIVE
  use Intrinsic, only: l_vmr
  implicit NONE
  private
  public :: ZATMOS_DERIV
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
  Subroutine ZATMOS_DERIV(ForwardModelConfig, FwdModelExtra, FwdModelIn, &
             frq_i, mid, delta, t_script, tau, ilo, ihi, k_atmos,Ier)
!
    type(forwardModelConfig_T), intent(in) :: forwardModelConfig
    type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra

    integer(i4), intent(in) :: mid, ILO, IHI, frq_i

    real(r8), intent(in) :: DELTA(:,:,:,:)
    real(r8), intent(in) :: T_SCRIPT(:), TAU(:)

    integer(i4), intent(out) :: Ier

    Type(path_derivative), intent(in out) :: k_atmos(:)
!
! --- Local variables --------------------------------------------
!
    Real(r8) :: r
    Integer(i4) :: J, IZ, IP
    type (VectorValue_T), pointer :: f
!
    Ier = 0
    do j = 1, size(forwardModelConfig%molecules)
!
      IF (forwardModelConfig%moleculeDerivatives(j)) THEN

        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
     &       quantityType=l_vmr, molecule=forwardModelConfig%molecules(j))
!
        do ip = 1, f%template%noInstances

          do iz = 1, f%template%noSurfs
!
            Call get_drad_notder (delta(1:,iz,ip,j), t_script, tau, &
   &                              mid, ilo, ihi, r)
            k_atmos(j)%values(frq_i,iz,ip) = r
!
          end do
!
        end do
!
      ENDIF
!
    end do
!
    Return
  End Subroutine ZATMOS_DERIV
end module ZATMOS_DERIV_M
