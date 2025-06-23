! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Do_Calc_Sparse_m

  implicit NONE
  private
  public :: Get_Do_Calc_Sparse
  public :: Get_Do_Calc_Sparse_Coarse !, Get_Do_Calc_Sparse_GL

  interface Get_Do_Calc_Sparse
    module procedure Get_Do_Calc_Sparse_Coarse !, Get_Do_Calc_Sparse_GL
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  subroutine Get_Do_Calc_Sparse_Coarse ( Tan_Pt_C, Eta, Col, N_Inds, Inds )

    ! Compute Inds(1:n_inds,1:2), indices on the coarse path for Col,
    ! such that one layer-boundary Eta is nonzero, or a GL Eta between
    ! the layer-boundaries is nonzero.

    use GLNP, ONLY: Ng, Ngp1
    use Sparse_m, only: Sparse_t

    integer, intent(in) :: Tan_Pt_C       ! Tangent point on the coarse path
    type(sparse_t), intent(in) :: Eta     ! Interpolating coefficients
    integer, intent(in) :: Col            ! Column of Eta_FZP
    integer,             intent(out) :: N_Inds    ! How many Inds produced
    integer, contiguous, intent(out) :: Inds(:,:) ! Nx2, Inds(:,1) is the
                                          ! layer index and the index of the
                                          ! boundary nearest the tangent point;
                                          ! Inds(:,2) is the other boundary index.

    logical :: Fill ! Fill next elements of inds
    integer :: NC   ! Next element of Eta in column Col
    integer :: Test

    n_inds = 0
    nc = eta%cols(col)       ! Last entry in the column
    if ( nc == 0 ) return    ! Column is empty
    do
      nc = eta%e(nc)%nc
      test = eta%e(nc)%r
      fill = n_inds <= 0
      if ( .not. fill ) fill = inds(n_inds,1) /= test ! Avoid duplicates
      if ( fill ) then
        n_inds = n_inds + 1
        inds(n_inds,1) = eta%e(nc)%r
        inds(n_inds,2) = inds(n_inds,1) + merge(-1,1,test <= tan_pt_c) * ngp1
      end if
      if ( nc == eta%cols(col) ) exit ! Just processed the last one
    end do

  end subroutine Get_Do_Calc_Sparse_Coarse

!----------------------------------------------------------------------
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, not_used_here ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Get_Do_Calc_Sparse_m

! $Log$
! Revision 2.1  2018/11/19 21:51:48  vsnyder
! Initial commit, in case it's ever useful for polarized molecule derivatives
!
