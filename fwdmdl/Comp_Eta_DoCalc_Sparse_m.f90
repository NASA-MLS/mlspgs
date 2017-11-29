! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Eta_DoCalc_Sparse_m

  implicit NONE

  private

  public :: Comp_Eta_Docalc_Sparse
  public :: Comp_All_Eta_2D,  Comp_One_Eta_2D
  public :: Comp_All_Eta_QTM, Comp_One_Eta_QTM
  public :: Get_Eta_DoCalc_ZP, Get_One_Eta_DoCalc_ZP

  interface Comp_Eta_Docalc_Sparse
    module procedure Comp_All_Eta_2D,  Comp_One_Eta_2D
    module procedure Comp_All_Eta_QTM, Comp_One_Eta_QTM
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------  Comp_All_Eta_2D  -----

  subroutine Comp_All_Eta_2D ( Grids_f, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP, Skip )

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f       ! Quantity values
    integer, intent(in) :: Tan_Pt              ! To split Z_Path into two
                                               ! monotone halves
    real(rp), intent(in) :: Z_Path(:)          ! Zetas on the path
    class(sparse_eta_t), intent(inout) :: Eta_Z(:) ! from VMR's zeta to path.
                                               ! InOut so as not to reallocate
                                               ! Eta_Z%Eta
    real(rp), intent(in) :: Phi_Path(:)        ! Phis on the path
    type(sparse_eta_t), intent(out) :: Eta_Phi(:)  ! Created here
    class(sparse_eta_t), intent(out) :: Eta_zP(:)  ! Created here
    integer, intent(in), optional :: Skip      ! At tangent point, default 1

    integer :: Col1   ! Last column in Eta_ZP before species I
    integer :: I

    col1 = 0
    do i = 1, size(grids_f%mol)
      call comp_one_eta_2d ( grids_f, i, tan_pt, z_path, eta_z(i), phi_path, &
                           & eta_phi(i), eta_zP(i), skip, col1 )
      col1 = col1 + ( grids_f%l_z(i) - grids_f%l_z(i-1) ) * &
                  & ( grids_f%l_p(i) - grids_f%l_p(i-1) )
    end do

  end subroutine Comp_All_Eta_2D

  ! -------------------------------------------  Comp_All_Eta_QTM  -----

  subroutine Comp_All_Eta_QTM ( Grids_f, Z_Path, Eta_Z,  Eta_Path, Eta_zQ )

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f           ! Quantity values
    real(rp), intent(in) :: Z_Path(:)              ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z(:)   ! from VMR's zeta to path.
    class(sparse_eta_t), intent(in) :: Eta_Path(:) ! from Phi to path for
                                                   ! all species, because there's
                                                   ! only one QTM.
    class(sparse_eta_t), intent(out) :: Eta_zQ(:)  ! Created here

    integer :: I

    do i = 1, size(grids_f%mol)
      call comp_one_eta_QTM ( grids_f, i, z_path, eta_z(i), &
                            & eta_path(i), eta_zQ(i) )
    end do

  end subroutine Comp_All_Eta_QTM

  ! --------------------------------------------  Comp_One_Eta_2D  -----

  subroutine Comp_One_Eta_2D ( Grids_f, N, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP, Skip, Col1 )

    use intrinsic, only: Lit_Indices
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f         ! Quantity values
    integer, intent(in) :: N                     ! Which quantity
    integer, intent(in) :: Tan_Pt                ! To split Z_Path into two
                                                 ! monotone halves
    real(rp), intent(in) :: Z_Path(:)            ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z    ! from VMR's Zeta to path.
    real(rp), intent(in) :: Phi_Path(:)          ! Phis on the path, Radians
    class(sparse_eta_t), intent(out) :: Eta_Phi  ! from VMR's Phi to path.
    class(sparse_eta_t), intent(inout) :: Eta_zP ! Created here
    integer, intent(in), optional :: Skip        ! At tangent point, default 1
    integer, intent(in), optional :: Col1        ! Value to put in Eta_zp%col1

    integer :: P1, P2, Z1, Z2 ! Boundaries from Grids_f%l_[pz]
    integer :: MySkip
    integer :: N_Path

    mySkip = 1
    if ( present(skip) ) mySkip = skip

    n_path = size(z_path) ! Assumed == size(p_path)

    p1 = grids_f%l_p(n-1)+1
    p2 = grids_f%l_p(n)
    z1 = grids_f%l_z(n-1)+1
    z2 = grids_f%l_z(n)

    ! Create and compute Eta_Phi
    call eta_phi%eta_1d ( grids_f%phi_basis(p1:p2), phi_path, &
                         & what=lit_indices(grids_f%mol(n)), resize=.true. )
    eta_phi%lbnd = [ p1 ] ! Column bounds
    eta_phi%ubnd = [ p2 ] !   in phi_basis

    ! Create Eta_Z
    call eta_z%create ( n_path, z2-z1+1, 2*(z2-z1+1), &
                      & ubnd = [ z2 ], lbnd = [ z1 ], &
                      & what=lit_indices(grids_f%mol(n)) )
    ! Compute the two halves of Eta_Z on either side of the tangent point
    ! separately.  These parts are monotone.  This avoids sorting Z_Path.
    call eta_z%eta_1d ( grids_f%zet_basis(z1:z2), z_path, &
                      & row1=tan_pt, rowN=1, create=.false. )
    call eta_z%eta_1d ( grids_f%zet_basis(z1:z2), z_path, &
                       & row1=tan_pt+mySkip, rowN=n_path, create=.false. )

    ! Compute Eta_ZP
    call eta_zp%eta_nd ( eta_z, eta_phi, &
                       & what=lit_indices(grids_f%mol(n)), resize=.true. )
    if ( present(col1) ) eta_zp%col1 = col1

  end subroutine Comp_One_Eta_2D

  ! -------------------------------------------  Comp_One_Eta_QTM  -----

  subroutine Comp_One_Eta_QTM ( Grids_f, N, Z_Path, Eta_Z, Eta_Path, Eta_zQ )

    use Sparse_Eta_m, only: Sparse_Eta_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f        ! Quantity values
    integer, intent(in) :: N                    ! Which quantity
    real(rp), intent(in) :: Z_Path(:)           ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z   ! from VMR's zeta to path.
    class(sparse_eta_t), intent(in) :: Eta_Path ! from Phi to path.
    class(sparse_eta_t), intent(out) :: Eta_zQ  ! Created here

!     integer :: I1, I2                         ! Boundaries from Grids_f%l_*
! 
!     eta_z%n = size(z_path,1)
!     i1 = grids_f%l_z(n-1)+1
!     i2 = grids_f%l_z(n)
!     ! Get the Zeta interpolator list
!     call get_eta_list ( grids_f%zet_basis(i1:i2), z_path, eta_z, sorted=.false. )
!     eta_zQ%n = eta_z%n
!     call get_eta_list ( eta_z%eta(1:eta_z%n), eta_path, eta_zQ )

  end subroutine Comp_One_Eta_QTM

  ! ------------------------------------------  Get_Eta_DoCalc_ZP  -----

  subroutine Get_Eta_DoCalc_ZP ( Eta_ZP_Sparse, Eta_ZP, Do_Calc_ZP, &
                               & NZ_ZP, NNZ_ZP, One_Sps )

  ! Compute Eta_ZP and Do_Calc_ZP from Eta_ZP_Sparse.  When the forward model
  ! no longer needs these, this subroutine can be deleted.

    use MLSKinds, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Sparse_Eta_m, only: Sparse_Eta_t

    class(sparse_eta_t), intent(in) :: Eta_ZP_Sparse(:) ! Created by Comp_*_Eta_* above
    real(rp), intent(inout) :: Eta_ZP(:,:) ! Eta_z * Eta_phi for each state
                           ! vector element. Marked inout so that we can clear
                           ! only the nonzeros without making the rest of it
                           ! undefined.
    logical, intent(inout), optional :: Do_Calc_ZP(:,:) ! Do_Calc_ZP(p,v)
                           ! indicates whether there is a contribution for 
                           ! state vector element V at point P on the path.
    integer, intent(inout), optional, target :: NZ_ZP(:,:) ! Nonzeros in eta_ZP
    integer, intent(inout), optional, target :: NNZ_ZP(:)  ! Number of nonzeros in eta_ZP
    integer, intent(in), optional :: One_Sps ! Do extraction for only this species

    integer :: J
    integer :: S, S1, S2   ! Species index = index in Eta_ZP_Sparse

    if ( present(nz_zp) .and. present(nnz_zp) ) then
      if ( present(do_calc_zp) ) then
        do j = 1, size(nnz_zp)
          do_calc_zp(nz_zp(1:nnz_zp(j),j),j) = .false.
        end do
      end if
      do j = 1, size(nnz_zp)
        eta_zp(nz_zp(1:nnz_zp(j),j),j) = 0
        nnz_zp(j) = 0
      end do
    else
      eta_zp = 0
      if ( present(do_calc_zp) ) do_calc_zp = .false.
    end if

    if ( present(one_sps) ) then
      s1 = one_sps
      s2 = one_sps
    else
      s1 = 1
      s2 = size(eta_zp_sparse,1)
    end if
    ! Get one less than the starting column in Eta_zp for species S.
    ! Columns have zeta-major order.  You might be tempted to put
    ! Grids_f%l_v(s-1) into Eta_zp_sparse%col1, but that fails for columns
    ! after species that have frequency dependence if Eta_zp_sparse isn't
    ! really Eta_fzp_sparse.
    do s = s1, s2
      if ( .not. allocated(eta_zp_sparse(s)%rows) .or. &
         & .not. allocated(eta_zp_sparse(s)%cols) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Rows or Cols of eta_zp_sparse(s) not allocated in Get_Eta_DoCalc_zp' )
        return
      end if
      call get_one_eta_docalc_zp ( eta_zp_sparse(s), eta_zp, do_calc_zp, &
                                 & nz_zp, nnz_zp )
    end do ! s

  end subroutine Get_Eta_DoCalc_ZP

  subroutine Get_One_Eta_DoCalc_ZP ( Eta_ZP_Sparse, Eta_ZP, Do_Calc_ZP, &
                                   & NZ_ZP, NNZ_ZP )

  ! Compute Eta_ZP and Do_Calc_ZP from Eta_ZP_Sparse.  When the forward model
  ! no longer needs these, this subroutine can be deleted.

    use Get_Do_Calc_m, only: Clean_Out_Nonzeros
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    class(sparse_eta_t), intent(in) :: Eta_ZP_Sparse ! Created by Comp_*_Eta_* above
    real(rp), intent(inout) :: Eta_ZP(:,:) ! Eta_z * Eta_phi for each state
                           ! vector element. Marked inout so that we can clear
                           ! only the nonzeros without making the rest of it
                           ! undefined.
    logical, intent(inout), optional :: Do_Calc_ZP(:,:) ! Do_Calc_ZP(p,v)
                           ! indicates whether there is a contribution for 
                           ! state vector element V at point P on the path.
    integer, intent(inout), optional, target :: NZ_ZP(:,:) ! Nonzeros in eta_ZP
    integer, intent(inout), optional, target :: NNZ_ZP(:)  ! Number of nonzeros in eta_ZP
    integer :: I, J, K
    integer :: L           ! Last column in Eta_ZP for species S-1

    call clean_out_nonzeros ( eta_zp, do_calc_zp, nz_zp, nnz_zp )
    l = eta_zp_sparse%col1
    do j = 1, size(eta_zp_sparse%cols)
      k = eta_zp_sparse%cols(j)
      if ( k == 0 ) cycle ! Column is empty
      i = k
      do
        i = eta_zp_sparse%e(i)%nc ! Next row in this column
        eta_zp ( eta_zp_sparse%e(i)%r,j+l ) = eta_zp_sparse%e(i)%v
        if ( present(nz_zp) ) then
          nnz_zp(j+l) = nnz_zp(j+l) + 1
          nz_zp(nnz_zp(j+l),j+l) = eta_zp_sparse%e(i)%r
        end if
        if ( present(do_calc_zp) ) do_calc_zp ( eta_zp_sparse%e(i)%r,j+l ) = .true.
        if ( i == k ) exit        ! No more rows in this column
      end do
    end do

  end subroutine Get_One_Eta_DoCalc_ZP

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Comp_Eta_DoCalc_Sparse_m

! $Log$
! Revision 2.5  2017/11/29 00:40:39  vsnyder
! Add Get_Eta_DoCalc_ZP, Get_One_Eta_DoCalc_ZP temporarily to get Eta matrix
! while full forward model still needs it.  Add "skip # at tangent" argument.
! Compute column starting position for each species.  Add thumbnails.  Use
! type-bound forms of Sparse_t procedures.
!
! Revision 2.4  2017/11/01 19:00:30  vsnyder
! Add tangent point
!
! Revision 2.3  2017/03/11 00:51:18  vsnyder
! Use Grids_F instead of Beta_Group
!
! Revision 2.2  2017/01/14 01:57:09  vsnyder
! Eliminate polymorphic interpolators.  Add arguments to return 1D Etas.
! Assume Z_Path and Phi_Path are not sorted.  Use template%Phi as radians
! because Phi_Path is radians.
!
! Revision 2.1  2016/12/15 02:44:32  vsnyder
! Initial commit
!
