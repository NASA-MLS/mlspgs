!
module ATMOS_BASIS_M
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component, &
                                  DEG2RAD
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP
  use L2PCDim, only: MNP => max_no_phi, NLVL
  use MLSCommon, only: I4, R4, R8
  use D_LINTRP_M, only: LINTRP
  use GET_ATMOS_M, only: GET_ATMOS
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-------------------------------------------------------------------

SUBROUTINE atmos_basis(time_stamp,atmospheric,no_atmos,z_grid,   &
     &                 spsfunc,n_lvls,mr_f,f_basis,no_coeffs_f,  &
     &                 phi_tan,no_phi_f,phi_basis_f,InDir,ld,ier)

!  ===============================================================
!  Declaration of variables for sub-program: atmos_basis
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Character (LEN=*), INTENT(IN) :: time_stamp, InDir

Integer(i4), INTENT(IN) :: no_atmos, n_lvls, ld

Integer(i4), INTENT(OUT) :: no_coeffs_f(:), no_phi_f(:), ier

Real(r8), INTENT(IN) :: phi_tan, z_grid(:)

Real(r8), INTENT(OUT) :: f_basis(:,:), phi_basis_f(:,:), mr_f(:,:,:), &
                         spsfunc(:,:)

type (atmos_comp), INTENT(IN) :: atmospheric(*)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: no_pts, ih, sps_i, jp, j, k, no_phi

Real(r8) :: fdum(nlvl,mnp), pdum(nlvl), dphi, tri_base(mxco)

! Begin code:

  dphi = 1.5 * deg2rad

  DO sps_i = 1, no_atmos

    CALL get_atmos(atmospheric(sps_i)%spectag,time_stamp,fdum, &
        pdum,no_pts,no_phi,InDir,ld,ier)
    IF(ier /= 0) RETURN

! Get number of coefficients on the retrieval grid

    k = atmospheric(sps_i)%no_lin_values
    no_coeffs_f(sps_i) = k
    no_phi_f(sps_i) = no_phi

! Create the phi_basis for this specie:

    ih = (no_phi + 1) / 2
    DO jp = 1, no_phi
      phi_basis_f(jp,sps_i) = phi_tan + (jp - ih) * dphi
    END DO

! Create mr_f:

    DO j = 1, k
      tri_base(j) = atmospheric(sps_i)%basis_peaks(j)
    END DO

! Loop over the phi's

    DO jp = 1, no_phi

! From the pdum grid, interpolate coeff. onto user's inputed atmospheric
! pressure grid, creating mr_f:

      CALL Lintrp(pdum,tri_base,fdum(1:,jp),mr_f(1:,jp,sps_i),no_pts,k)

    END DO

!  *** TESTING CODE: Make all mixing ratio equal to the "center" Phi value
!  (Thus no gradients in mr_f in the "Phi" dimension..)

    DO jp = 1, no_phi
      IF(jp /= ih) THEN
        DO j = 1, k
          mr_f(j,jp,sps_i) = mr_f(j,ih,sps_i)
        END DO
      END IF
    END DO

!  *** END TESTING CODE

! Copy user pressure grid into f_basis grid:

  DO j = 1, k
    f_basis(j,sps_i) = tri_base(j)
  END DO

! Interpolate atmospheric function (spsfunc) onto the preselected grid

    CALL Lintrp(tri_base,z_grid,mr_f(1:,ih,sps_i),spsfunc(1:,sps_i),k,n_lvls)

  END DO

  RETURN
END SUBROUTINE atmos_basis
end module ATMOS_BASIS_M
! $Log$
! Revision 1.1  2000/06/21 21:56:09  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
