! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module testfield_m

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!-----------------------------------------------------------------------

contains 

subroutine testfield ( h, ct, stcp, stsp)
    use MLSCommon, only: Rk => Rp
    real(rk), intent(inout) ::    h(:)      ! mag field magnitude
    real(rk), intent(inout) ::    CT(:)     ! Cos(Theta)
    real(rk), intent(inout) ::    STCP(:)   ! Sin(Theta) Cos(Phi)
    real(rk), intent(inout) ::    STSP(:)   ! Sin(Theta) Sin(Phi)

    h=0
    ct=0
    stcp=0.707
    stsp=0.707
end subroutine testfield
end module testfield_m
! $Log$
! Revision 1.1.2.1  2003/03/19 00:17:53  michael
! Procedure to cobble dummy magnetic field path into code without recompiling
! all of Full_Foreward_Model
!
