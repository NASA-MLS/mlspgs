subroutine GET_3D_VIEW ( Input, I1, I2, I3, Output )
  ! The purpose of this routine is to "cheat" to get a rank-3 view of a
  ! (piece of a) rank-1 vector.  It's declared here with "input" having
  ! rank 3.  It is intended that at its point of use, it will be declared
  ! (by way of an interface body) so that "input" is rank-1 assumed size
  ! (i.e. with its dimension = *).  Therefore, this subroutine cannot be
  ! a module procedure -- the compiler would catch us at our lying.
  ! Hopefully, the Fortran standard will someday allow us to do this
  ! harmless and useful operation directly.
  use MLSCommon, only: R8
  integer, intent(in) :: I1, I2, I3
  real(r8), intent(in), target :: Input(I1, I2, I3)
  real(r8), pointer :: Output(:,:,:)
  output => input
end subroutine GET_3D_VIEW
