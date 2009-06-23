! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GET_DRAD_NOTDER_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: GET_DRAD_NOTDER
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
! This subroutine computes the scalarized condensed radiative transfer
! derivatives, without Temperature derivatives
!
  Subroutine GET_DRAD_NOTDER (DDER, T_SCRIPT, TAU, mid, ILO, IHI, RDRAD)
!
    Integer(i4), intent(in) :: mid, ILO, IHI
    Real(r8), intent(in) :: DDER(:), T_SCRIPT(:), TAU(:)
    Real(r8), intent(out) :: RDRAD

    Integer(i4) :: h_i, IndxL
    Real(r8) :: drad, q, w
!
    w = 0.0_r8
    drad = 0.0_r8
!
    h_i = 2
    do while (h_i <= ilo)
      w = w - dder(h_i-1)
      q = t_script(h_i) * w
      drad = drad + q * tau(h_i)
      h_i = h_i + 1
    end do
!
    if (ilo > ihi) then
      rdrad = drad
      Return
    endif
!
    IndxL = mid + 1
    q = t_script(IndxL) * w
    drad = drad + q * tau(IndxL)
!
    h_i = IndxL + 1
    do while (h_i <= ihi)
      w = w - dder(h_i-2)
      q = t_script(h_i) * w
      drad = drad + q * tau(h_i)
      h_i = h_i + 1
    end do
!
    rdrad = drad
!
    Return
  End Subroutine GET_DRAD_NOTDER
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module GET_DRAD_NOTDER_M
! $Log$
! Revision 2.2  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.2  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.1  2001/03/31 23:18:08  zvi
! Initial release - get_drad w/o Temp. derivatives
!
! Revision 1.1  2001/03/31 18:12:05  Z.Shippony
! Initial conversion to Fortran 90
