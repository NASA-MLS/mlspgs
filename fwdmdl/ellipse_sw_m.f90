! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ELLIPSE_SW_M
  use MLSCommon, only: R8
  use ELLIPSE_M, only: ELLIPSE
  Implicit NONE
  Private
  Public :: PHI_TO_H_S, H_TO_S_PHI, S_TO_H_PHI

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
!----------------------------------------------------  Phi_to_H_S  -----
!  Given Phi, get H and S
!  ** Note: This routine is using The Equivalent Circle concept

  Subroutine Phi_to_H_S ( Elvar, Phi, H, S )

    Real(r8), Intent(in) :: Phi
    Real(r8), Intent(out) :: H, S

    Type(ellipse), intent(in out) :: Elvar

    Real(r8) :: Delphi, Q, Rt, V

    delphi = Phi - elvar%Phi_tan
    v = Cos(delphi)

    if ( .not. elvar%EarthX ) then
      rt = elvar%ht + elvar%RoC
      h = rt / v - elvar%RoC
      S = rt * Tan(delphi)
    else
      if ( elvar%ps > 0.0 ) &
        &  v = Cos(Phi-(2.0*elvar%Phi_s-elvar%Phi_tan))
      h = elvar%RoC * (elvar%Rr / v - 1.0)
      q = Sin(Phi-elvar%Phi_s) / v
      S = elvar%RoC * abs(q)
    end if

    Return
  End Subroutine Phi_to_H_S

! ---------------------------------------------------  H_to_S_Phi  -----
!  Given H, get S and Phi
!  ** Note: This routine is using The Equivalent Circle concept

  Subroutine H_to_S_Phi ( Elvar, H, S, Phi )

    Real(r8), intent(in) :: H
    Real(r8), intent(out) :: S, Phi

    Type(ellipse), intent(in out) :: Elvar

    Real(r8) :: Q, R, Rt, V

    r = h + elvar%RoC

    if ( .not. elvar%EarthX ) then
      S = 0.0
      Phi = elvar%Phi_tan
      rt = elvar%ht + elvar%RoC
      q = r * r - rt * rt
      if ( q >= 1.0e-6_r8 ) S = Sqrt(q)
      v = rt / r
      if ( abs(v) > 1.0 ) v = Sign(1.0_r8,v)
      Phi = elvar%Phi_tan + elvar%ps * Acos(v)
    else
      v = (elvar%RoC / r) * elvar%Rr   ! Rr = Cos(Phi_tan-Phi_s)=(ht+RoC)/RoC
      if ( abs(v) > 1.0 ) v = Sign(1.0_r8,v)
      if ( elvar%ps < 0.0 ) then
        Phi = elvar%Phi_tan - Acos(v)
      else
        Phi = 2.0 * elvar%Phi_s - elvar%Phi_tan + Acos(v)
      end if
      q = Sin(Phi-elvar%Phi_s) / v
      S = elvar%RoC * abs(q)
    end if

    Return
  End Subroutine H_to_S_Phi

!----------------------------------------------------  S_to_H_Phi  -----
!  Given S, get H and Phi
!  ** Note: This routine is using The Equivalent Circle concept

  Subroutine S_to_H_Phi ( Elvar, S, H, Phi )

    Real(r8), intent(in) :: S
    Real(r8), intent(out) :: Phi, H

    Type(ellipse), intent(in out) :: ELVAR

    Real(r8) :: Q, R, Rt, V

    if ( .not. elvar%EarthX ) then
      rt = elvar%ht + elvar%RoC
      r = Sqrt(S*S + rt*rt)
      h = r - elvar%RoC
      v = rt / r
      if ( abs(v) > 1.0 ) v = Sign(1.0_r8,v)
      Phi = elvar%Phi_tan + elvar%ps * Acos(v)
    else
      q = elvar%ps * S / elvar%RoC
      if ( elvar%ps < 0.0 ) then
        Phi = Atan2(q*elvar%cpt+elvar%sps,elvar%cps-q*elvar%spt)
        v = Cos(Phi - elvar%Phi_tan)
      else
        Phi = Atan2(q*elvar%cpts+elvar%sps,elvar%cps-q*elvar%spts)
        v = Cos(Phi- 2.0 * elvar%Phi_s + elvar%Phi_tan)
      end if
      h = elvar%RoC * (elvar%Rr / v - 1.0)
    end if

    Return
  End Subroutine S_to_H_Phi
!---------------------------------------------------------------------------
end module ELLIPSE_SW_M

! $Log$
! Revision 1.5  2001/04/19 22:54:46  vsnyder
! Use generic ACOS instead of specific DACOS
!
! Revision 1.4  2001/04/19 06:48:13  zvi
! Fixing memory leaks..
!
! Revision 1.3  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.1  2000/05/04 18:12:05  Z.Shippony
! Initial conversion to Fortran 90
!
