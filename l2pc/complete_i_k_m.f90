module COMPLETE_I_K_M
  use L2PCDim, only: MNP => max_no_phi, NLVL, NSPS
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PC_PFA_STRUCTURES
  implicit NONE
  private
  public :: COMPLETE_I_K
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Subroutine COMPLETE_I_K ( J, NO_HTS, NO_GEOPHYS, I_RAW, K_STAR_ATMOS, &
 &                          K_STAR_GEOPHYS )
!  Copy the 'j' values onto the 'j+1,j+2,...No_Hts' locations in all
!  arrays and matrices ..
    Integer(i4), intent(in) :: J
    Integer(i4), intent(in) :: NO_HTS
    Integer(i4), intent(in) :: NO_GEOPHYS
    Real(r8), intent(inout) :: I_RAW(*)
    Real(r4), intent(inout) :: K_STAR_ATMOS(Nlvl,mxco,mnp,*)
    Real(r4), intent(inout) :: K_STAR_GEOPHYS(Nlvl,mxco,mnp,*)
    Integer(i4) :: I,K,N                ! Subscripts and loop inductors
!  Complete the radiances last location (No_Hts)
    i_raw(j+1:No_Hts) = i_raw(j)
!  Complete the geophysical FIRST derivative matrices last location
    do k = 1, no_geophys
      do i = 1, mnp
        do n = 1, mxco
          k_star_geophys(j+1:No_Hts,n,i,k) = k_star_geophys(j,n,i,k)
        end do
      end do
    end do
!  Complete the atmospheric FIRST derivative matrices last location
    do k = 1, Nsps
      do i = 1, mnp
        do n = 1, mxco
          k_star_atmos(j+1:No_Hts,n,i,k) = k_star_atmos(j,n,i,k)
        end do
      end do
    end do
    Return
  End Subroutine COMPLETE_I_K
end module COMPLETE_I_K_M
! $Log$
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
