! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

module Piq_int_m

  implicit none

  private
  public :: Piq_int
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character ( len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  contains
!---------------------------------------------------------------------------

  subroutine Piq_int ( z_grid, t_basis, z_ref, piq )

! Compute the piq (sans mass) used in the L2PC
! hydrostatic function
! This assumes the mean molecular mass is invariant
! between z_grid levels
! uses L2PC style triangular temperature representation basis
! argument z_ref allows a reference pressure that is /= z_grid(1)
! however adding this feature will cause program to run twice as slow

    use MLSCommon, only: RP, IP

    real(rp), intent(in) :: z_grid(:), t_basis(:), z_ref
    real(rp), intent(out) :: piq(:,:)

! inside code variables

    real(rp) :: a,c,aa,cc
    real(rp), dimension(1:size(z_grid)) :: b, d, bb, dd
    integer(ip) :: n_coeffs, i, ind(1)

! begin code
! Establish dimensions

    n_coeffs = size(t_basis)

! locate z_ref relative to t_basis
! I don't know if this is the fastest way to do this but this is a
! method derived from idl's ind_scl routine

    ind = min(max(pack((/(i,i=0,n_coeffs)/), &
                       (/minval((/z_ref,(t_basis(i),i=1,n_coeffs)/))-1.0_rp, &
                          (t_basis(i),i=1,n_coeffs)/) <= z_ref .and. &
                          z_ref < (/(t_basis(i),i=1,n_coeffs), &
                       maxval((/z_ref,(t_basis(i),i=1,n_coeffs)/))+1.0_rp/)), &
                  2), &
               n_coeffs-2)

! NOTE it seems like MIN/MAX need to be done in pairs when comparing
! against arrays and scalers.
! for all coeffients below ind use
! initial coefficient

    a = max(t_basis(1),z_ref)
    b = max(min(max(z_grid,t_basis(1)),t_basis(2)),z_ref)
    c = min(z_ref,t_basis(2))
    d = min(min(max(z_grid,t_basis(1)),z_ref),t_basis(2))
    piq(:,1) = max(min(z_grid,t_basis(1)),z_ref) - z_ref &
             + min(min(z_grid,t_basis(1)),z_ref) - min(t_basis(1),z_ref) &
             + ((t_basis(2)-0.5_rp*(a+b))*(b-a) &
             -  (t_basis(2)-0.5_rp*(c+d))*(c-d))/(t_basis(2)-t_basis(1))

! these are wholly negative contributions only

    do i = 2,ind(1) - 1
      b = min(max(z_grid,t_basis(i-1)),t_basis(i))
      d = min(max(z_grid,t_basis(i)),t_basis(i+1))
      piq(:,i) = (t_basis(i) - b)*(0.5_rp*(t_basis(i) - b) &
               / (t_basis(i) - t_basis(i-1)) - 1.0_rp) &
               - 0.5_rp*(t_basis(i+1) - d)**2 / (t_basis(i+1)-t_basis(i))
    end do

! coefficients where z_ref is amongst t_basis

    do i = ind(1),ind(1)+1
      a = max(t_basis(i-1),z_ref)
      b = max(min(max(z_grid,t_basis(i-1)),t_basis(i)),z_ref)
      c = min(z_ref,t_basis(i))
      d = min(min(max(z_grid,t_basis(i-1)),z_ref),t_basis(i))
      aa = max(t_basis(i),z_ref)
      bb = max(min(max(z_grid,t_basis(i)),t_basis(i+1)),z_ref)
      cc = min(z_ref,t_basis(i+1))
      dd = min(min(max(z_grid,t_basis(i)),z_ref),t_basis(i+1))
      piq(:,i) = ((0.5_rp*(a+b)-t_basis(i-1))*(b-a) &
               -  (0.5_rp*(c+d)-t_basis(i-1))*(c-d)) &
               / (t_basis(i)-t_basis(i-1)) &
               + ((t_basis(i+1)-0.5_rp*(aa+bb))*(bb-aa) &
               -  (t_basis(i+1)-0.5_rp*(cc+dd))*(cc-dd)) &
               / (t_basis(i+1)-t_basis(i))
    end do

! for all coeffients above ind use

    do i = ind(1)+2,n_coeffs-1
      b = min(max(z_grid,t_basis(i-1)),t_basis(i))
      d = min(max(z_grid,t_basis(i)),t_basis(i+1))
      piq(:,i) = 0.5_rp*(b - t_basis(i-1))**2 / (t_basis(i)-t_basis(i-1)) &
               + (d - t_basis(i))*(1.0 - 0.5_rp*(d - t_basis(i)) &
               / (t_basis(i+1)-t_basis(i)))
    end do

! upper coefficient

    a = max(t_basis(n_coeffs-1),z_ref)
    b = max(min(max(z_grid,t_basis(n_coeffs-1)),t_basis(n_coeffs)),z_ref)
    c = min(z_ref,t_basis(n_coeffs))
    d = min(min(max(z_grid,t_basis(n_coeffs-1)),z_ref),t_basis(n_coeffs))
    piq(:,n_coeffs) = ((0.5_rp*(a+b)-t_basis(n_coeffs-1))*(b-a) &
                    -  (0.5_rp*(c+d)-t_basis(n_coeffs-1))*(c-d)) &
                    / (t_basis(n_coeffs)-t_basis(n_coeffs-1)) &
                    + max(max(z_grid,t_basis(n_coeffs)),z_ref) &
                    - max(t_basis(n_coeffs),z_ref) &
                    - min(z_ref,max(z_grid,t_basis(n_coeffs))) &
                    + min(t_basis(n_coeffs),z_ref)

  end subroutine Piq_int

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Piq_int_m
!---------------------------------------------------
! $Log$
! Revision 2.1  2002/09/25 20:35:30  vsnyder
! Move USE from module scope to procedure scope.  Change allocatable arrays
! to automatic arrays.  Insert the copyright notice.  Cosmetic changes.
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1.2.3  2001/09/13 22:51:23  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:52  zvi
! Added CVS stuff
!
