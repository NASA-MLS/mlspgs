!
module ELLIPSE_SW_M
  use MLSCommon, only: R8
  use ELLIPSE_M, only: ELLIPSE
  Implicit NONE
  Private
  Public PHI_TO_H_S, H_TO_S_PHI, S_TO_H_PHI
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
!  Given Phi, get H and S
!  ** Note: This routine is using The Equivalent Circel concept
!
  Subroutine Phi_to_H_S(elvar,Phi,h,S)
!
    Real(r8), Intent(IN) :: Phi
    Real(r8), Intent(OUT) :: H, S
!   
    Type(ELLIPSE), intent(in out) :: elvar
!
!
    Real(r8) :: rt,delphi,q,v
! 
    delphi = Phi - elvar%Phi_tan
    v = Cos(delphi)
!
    if(.not.elvar%EarthX) then
      rt = elvar%ht + elvar%RoC
      h = rt / v - elvar%RoC
      S = rt * Tan(delphi)
    else
      if(elvar%ps.gt.0.0) &
     &     v = Cos(Phi-(2.0d0*elvar%Phi_s-elvar%Phi_tan))
      h = elvar%RoC * (elvar%Rr / v - 1.0d0)
      q = Sin(Phi-elvar%Phi_s) / v
      S = elvar%RoC * abs(q)
    end if
!
    Return
  End Subroutine Phi_to_H_S
!
!------------------------------------------------------------------------
!  Given H, get S and Phi
!  ** Note: This routine is using The Equivalent Circel concept
!
      Subroutine H_to_S_Phi(elvar,h,S,Phi)
!
      Real(r8), intent(IN) :: h
      Real(r8), intent(OUT) :: S, Phi
!
      Type(ELLIPSE), intent(in out) :: elvar

      Real(r8) :: q,v,r,rt
!
      r = h + elvar%RoC
!
      if(.not.elvar%EarthX) then
        S = 0.0d0
        Phi = elvar%Phi_tan
        rt = elvar%ht + elvar%RoC
        q = r * r - rt * rt
        if(q.ge.1.0d-6) S = Sqrt(q)
        v = rt / r
        if(abs(v).gt.1.0d0) v = Sign(1.0_r8,v)
        Phi = elvar%Phi_tan + elvar%ps * Acos(v)
      else
        v = (elvar%RoC / r) * elvar%Rr   ! Rr = Cos(Phi_tan-Phi_s)=(ht+RoC)/RoC
        if(abs(v).gt.1.0d0) v = Sign(1.0_r8,v)
        if(elvar%ps.lt.0.0) then
          Phi = elvar%Phi_tan - Acos(v)
        else
          Phi = 2.0d0 * elvar%Phi_s - elvar%Phi_tan + Acos(v)
        end if
        q = Sin(Phi-elvar%Phi_s) / v
        S = elvar%RoC * abs(q)
      end if
!
      Return
      End Subroutine H_to_S_Phi
!
!------------------------------------------------------------------------
!  Given S, get H and Phi
!  ** Note: This routine is using The Equivalent Circel concept
!
      Subroutine S_to_H_Phi(elvar,S,h,Phi)
!
      Real(r8), intent(IN) :: S
      Real(r8), intent(OUT) :: Phi,H
!
      Type(ELLIPSE), intent(in out) :: elvar

      Real(r8) :: rt,r,q,v
!
      if(.not.elvar%EarthX) then
        rt = elvar%ht + elvar%RoC
        r = Sqrt(S*S + rt*rt)
        h = r - elvar%RoC
        v = rt / r
        if(abs(v).gt.1.0d0) v = Sign(1.0_r8,v)
        Phi = elvar%Phi_tan + elvar%ps * Acos(v)
      else
        q = elvar%ps * S / elvar%RoC
        if(elvar%ps.lt.0.0) then
          Phi = DAtan2(q*elvar%cpt+elvar%sps,elvar%cps-q*elvar%spt)
          v = Cos(Phi - elvar%Phi_tan)
        else
          Phi = DAtan2(q*elvar%cpts+elvar%sps,elvar%cps-q*elvar%spts)
          v = Cos(Phi- 2.0d0 * elvar%Phi_s + elvar%Phi_tan)
        end if
        h = elvar%RoC * (elvar%Rr / v - 1.0d0)
      end if
!
      Return
      End Subroutine S_to_H_Phi
!---------------------------------------------------------------------------
end module ELLIPSE_SW_M
! $Log$
! Revision 1.4  2001/04/19 06:48:13  zvi
! Fixing memory leaks..
!
! Revision 1.3  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.1  2000/05/04 18:12:05  Z.Shippony
! Initial conversion to Fortran 90
!
