module GH06
  use MLSCommon, only: R8
  private
  public :: GX, GW

  integer, private, parameter :: N = 3

  real(r8), parameter :: GX(n) = &
  & (/ 4.36077411927616508271d-1, 1.33584907401369694957d+0, &
  &    2.35060497367449222280d+0 /)
!
  real(r8), parameter :: GW(n) = &
  & (/ 7.24629595224392524608d-1, 1.57067320322856644036d-1, &
  &    4.53000990550884564224d-3 /)

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

end module GH06

! Log: GH06,v $
