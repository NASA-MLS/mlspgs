module GET_BETA_PATH_M
  use MLSCommon, only: I4, R8
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC
  use PATH_ENTITIES_M, only: PATH_VECTOR
  use CREATE_BETA_M, only: CREATE_BETA
  implicit NONE
  private
  public :: get_beta_path
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains

  Subroutine get_beta_path (Spectag, Freq, no_ele, z_path, t_path, &
                            mdb_hdr, mdb_rec, values, t_power, dbeta_dw, &
                            dbeta_dn, dbeta_dnu, Ier)

    integer(i4), intent(in) :: SPECTAG, no_ele

    real(r8), intent(in) :: Freq

    Type(path_vector), intent(in) :: z_path, t_path

    Type (eos_mdb_hdr), intent(in) :: MDB_HDR
    Type (eos_mdb_rec), intent(in) :: MDB_REC(*)

    integer(i4), intent(out) :: Ier
    real(r8), intent(out) :: values(*), t_power(*), dbeta_dw(*), &
                             dbeta_dn(*), dbeta_dnu(*)


! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: h_i
    Real(r8) :: pn, pnu, pw, q, t, t_p, z

! -----     Executable statements     ----------------------------------
!
!
    values(:no_ele) = 0.0
    t_power(:no_ele) = 0.0
    dbeta_dw(:no_ele) = 0.0
    dbeta_dn(:no_ele) = 0.0
    dbeta_dnu(:no_ele) = 0.0

    do h_i = 1, no_ele
      z = z_path%values(h_i)
      if(z <= -4.5) CYCLE
      t = t_path%values(h_i)
      Call Create_beta(Spectag, z, t, Freq, mdb_hdr, mdb_rec, &
   &                   q, t_p, pw, pn, pnu, Ier)
      if(Ier /= 0) Return
      values(h_i) = q
      t_power(h_i) = t_p
      dbeta_dw(h_i) = pw
      dbeta_dn(h_i) = pn
      dbeta_dnu(h_i) = pnu
    end do
!
    Return
!
  end subroutine get_beta_path
!
end module GET_BETA_PATH_M
! $Log$
! Revision 1.1  2000/06/21 21:56:10  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
