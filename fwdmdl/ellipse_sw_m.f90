!
module ELLIPSE_SW_M
  use MLSCommon, only: R8
  use ELLIPSE, only: CPT, SPT, CPS, SPS, CPTS, SPTS, HT, RR, PHI_TAN, &
                     PHI_S, PS, ROC, EARTHX
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
      Subroutine Phi_to_H_S(Phi,h,S)
!
      Real(r8), Intent(IN) :: Phi
      Real(r8), Intent(OUT) :: H, S
!
      Real(r8) :: rt,delphi,q,v
!
      delphi = Phi - Phi_tan
      v = Cos(delphi)
!
      if(.not.EarthX) then
        rt = ht + RoC
        h = rt / v - RoC
        S = rt * Tan(delphi)
      else
        if(ps.gt.0.0) v = Cos(Phi-(2.0d0*Phi_s-Phi_tan))
        h = RoC * (Rr / v - 1.0d0)
        q = Sin(Phi-Phi_s) / v
        S = RoC * abs(q)
      endif
!
      Return
      End Subroutine Phi_to_H_S
!
!------------------------------------------------------------------------
!  Given H, get S and Phi
!  ** Note: This routine is using The Equivalent Circel concept
!
      Subroutine H_to_S_Phi(h,S,Phi)
!
      Real(r8), intent(IN) :: h
      Real(r8), intent(OUT) :: S, Phi

      Real(r8) :: q,v,r,rt
!
      r = h + RoC
!
      if(.not.EarthX) then
        S = 0.0d0
        Phi = Phi_tan
        rt = ht + RoC
        q = r * r - rt * rt
        if(q.ge.1.0d-6) S = Sqrt(q)
        v = rt / r
        if(abs(v).lt.1.0d0) Phi = Phi_tan + ps * DAcos(v)
      else
        v = (RoC / r) * Rr        ! Rr = Cos(Phi_tan-Phi_s)=(ht+RoC)/RoC
        if(ps.lt.0.0) then
          Phi = Phi_tan - DAcos(v)
        else
          Phi = 2.0d0 * Phi_s - Phi_tan + DAcos(v)
        endif
        q = Sin(Phi-Phi_s) / v
        S = RoC * abs(q)
      endif
!
      Return
      End Subroutine H_to_S_Phi
!
!------------------------------------------------------------------------
!  Given S, get H and Phi
!  ** Note: This routine is using The Equivalent Circel concept
!
      Subroutine S_to_H_Phi(S,h,Phi)
!
      Real(r8), intent(IN) :: S
      Real(r8), intent(OUT) :: Phi,H

      Real(r8) :: rt,r,q,v
!
      if(.not.EarthX) then
        rt = ht + RoC
        r = Sqrt(S*S + rt*rt)
        h = r - RoC
        v = rt / r
        Phi = Phi_tan + ps * DAcos(v)
      else
        q = ps * S / RoC
        if(ps.lt.0.0) then
          Phi = DAtan2(q*cpt+sps,cps-q*spt)
          v = Cos(Phi - Phi_tan)
        else
          Phi = DAtan2(q*cpts+sps,cps-q*spts)
          v = Cos(Phi- 2.0d0 * Phi_s + Phi_tan)
        endif
        h = RoC * (Rr / v - 1.0d0)
      endif
!
      Return
      End Subroutine S_to_H_Phi
!---------------------------------------------------------------------------
end module ELLIPSE_SW_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  Z.Shippony
! Initial conversion to Fortran 90
!
