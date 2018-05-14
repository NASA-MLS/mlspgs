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
  public :: Comp_One_Eta_Z
  public :: Get_Eta_DoCalc_FZP
  ! Be careful with Get_One_Eta_DoCalc_ZP to send only the columns of
  ! Eta, Do_Calc, NZ and NNZ that are gemane to the species.
  public :: Get_One_Eta_DoCalc_FZP

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
                             & Phi_Path, Eta_Phi, Eta_zp, Eta_fzp, Skip )

  ! Compute interpolation coefficients for Zeta x Phi for all molecules

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
    type(sparse_eta_t), intent(inout) :: Eta_Phi(:)  ! Created here
    class(sparse_eta_t), intent(inout) :: Eta_zp(:)  ! Created here
    class(sparse_eta_t), intent(inout) :: Eta_fzp(:) ! Created here
    integer, intent(in), optional :: Skip      ! At tangent point, default 1

    integer :: I

    do i = 1, size(grids_f%mol)
      if ( grids_f%l_f(i-1)+1 == grids_f%l_f(i) ) then ! Not frequency dependent
        call comp_one_eta_2d ( grids_f, i, tan_pt, z_path, eta_z(i), &
                             & phi_path, eta_phi(i), eta_fzp(i), skip )
      else ! Frequency dependent, fill Eta_ZP so we can calculate Eta_FZP
           ! when we have frequency
        call comp_one_eta_2d ( grids_f, i, tan_pt, z_path, eta_z(i), &
                             & phi_path, eta_phi(i), eta_zp(i), skip )
      end if
    end do

  end subroutine Comp_All_Eta_2D

  ! -------------------------------------------  Comp_All_Eta_QTM  -----

  subroutine Comp_All_Eta_QTM ( Grids_f, Tan_Pt, Z_Path, Eta_Z, Eta_Path, &
                              & Eta_zQ, Skip )

  ! Compute interpolation coefficients for Zeta x QTM for all molecules

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f           ! Quantity values
    integer, intent(in) :: Tan_Pt                  ! To split Z_Path into two
                                                   ! monotone halves
    real(rp), intent(in) :: Z_Path(:)              ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z(:)   ! from VMR's zeta to path.
    class(sparse_eta_t), intent(in) :: Eta_Path(:) ! from QTM to path for
                                                   ! all species, because
                                                   ! there's only one QTM.
    class(sparse_eta_t), intent(out) :: Eta_zQ(:)  ! Created here
    integer, intent(in), optional :: Skip          ! At tangent point, default 1

    integer :: I

    ! Compute Eta_Path for all molecules because there's only one QTM
    ! for everything

    do i = 1, size(grids_f%mol)
      ! Compute Eta_ZQ(i) = Eta_Z(i) * Eta_Path(i)
      call comp_one_eta_QTM ( grids_f, i, tan_pt, z_path, eta_z(i), &
                            & eta_path(i), eta_zQ(i), skip )
    end do

  end subroutine Comp_All_Eta_QTM

  ! --------------------------------------------  Comp_One_Eta_2D  -----

  subroutine Comp_One_Eta_2D ( Grids_f, N, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP, Skip )

  ! Compute interpolation coefficients for Zeta x Phi for one molecule

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

    integer :: P1, P2 ! Boundaries from Grids_f%l_p
    integer :: What

    what = grids_f%qtyStuff(n)%qty%template%name

    p1 = grids_f%l_p(n-1)+1
    p2 = grids_f%l_p(n)

    ! Create and compute Eta_Phi
    call eta_phi%eta_1d ( grids_f%phi_basis(p1:p2), phi_path, &
                        & what=what, resize=.true. )
    eta_phi%lbnd = [ p1 ] ! Column bounds
    eta_phi%ubnd = [ p2 ] !   in phi_basis

    ! Create Eta_Z
    call comp_one_eta_z ( grids_f, n, tan_pt, z_path, eta_z, skip )

    ! Compute Eta_ZP
    call eta_zp%eta_nd ( eta_z, eta_phi, what=what, resize=.true. )

  end subroutine Comp_One_Eta_2D

  ! -------------------------------------------  Comp_One_Eta_QTM  -----

  subroutine Comp_One_Eta_QTM ( Grids_f, N, Tan_Pt, Z_Path, Eta_Z, Eta_Path, &
                              & Eta_zQ, Skip )

  ! Compute interpolation coefficients for Zeta x QTM for one molecule

    use Sparse_Eta_m, only: Sparse_Eta_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f        ! Quantity values
    integer, intent(in) :: N                    ! Which quantity
    integer, intent(in) :: Tan_Pt               ! To split Z_Path into two
                                                ! monotone halves
    real(rp), intent(in) :: Z_Path(:)           ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z   ! from VMR's zeta to path.
    class(sparse_eta_t), intent(in) :: Eta_Path ! from QTM to path.  Not
                                                ! computed here because there's
                                                ! only one QTM for everything
    class(sparse_eta_t), intent(out) :: Eta_zQ  ! Created here
    integer, intent(in), optional :: Skip       ! At tangent point, default 1

    ! Create Eta_Z
    call comp_one_eta_z ( grids_f, n, tan_pt, z_path, eta_z, skip )

    ! Compute Eta_ZQ
!     call eta_zp%eta_nd ( eta_z, eta_path, what=what, resize=.true. )

  end subroutine Comp_One_Eta_QTM

  ! ---------------------------------------------  Comp_One_Eta_Z  -----

  subroutine Comp_One_Eta_Z ( Grids_f, N, Tan_Pt, Z_Path, Eta_Z, &
                            & Skip )

  ! Compute interpolation coefficients for Zeta for one molecule

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f         ! Quantity values
    integer, intent(in) :: N                     ! Which quantity
    integer, intent(in) :: Tan_Pt                ! To split Z_Path into two
                                                 ! monotone halves
    real(rp), intent(in) :: Z_Path(:)            ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z    ! from VMR's Zeta to path.
    integer, intent(in), optional :: Skip        ! At tangent point, default 1

    integer :: Z1, Z2 ! Boundaries from Grids_f%l_z
    integer :: MySkip
    integer :: N_Path
    integer :: What

    what = grids_f%qtyStuff(n)%qty%template%name

    mySkip = 1
    if ( present(skip) ) mySkip = skip

    n_path = size(z_path) ! Assumed == size(p_path)

    z1 = grids_f%l_z(n-1)+1
    z2 = grids_f%l_z(n)

    ! Create Eta_Z
    call eta_z%create ( n_path, z2-z1+1, 2*(z2-z1+1), &
                      & ubnd = [ z2 ], lbnd = [ z1 ], what=what )
    ! Compute the two halves of Eta_Z on either side of the tangent point
    ! separately.  These parts are monotone.  This avoids sorting Z_Path.
    call eta_z%eta_1d ( grids_f%zet_basis(z1:z2), z_path, &
                      & row1=tan_pt, rowN=1, create=.false. )
    call eta_z%eta_1d ( grids_f%zet_basis(z1:z2), z_path, &
                      & row1=tan_pt+mySkip, rowN=n_path, create=.false. )

  end subroutine Comp_One_Eta_Z

  ! -----------------------------------------  Get_Eta_DoCalc_FZP  -----

  subroutine Get_Eta_DoCalc_FZP ( Eta_FZP_Sparse, Eta_FZP, Grids_f, &
                                & Do_Calc_FZP, NZ_FZP, NNZ_FZP, One_Sps, &
                                & Eta_ZP_Sparse, Spread )

  ! Compute Eta_FZP and Do_Calc_FZP from Eta_FZP_Sparse and Eta_ZP_Sparse.
  ! Where Grids_f%deriv_flags is false, don't set Do_Calc_FZP 
  ! When the forward model no longer needs these, this subroutine can be
  ! deleted.

    use Comp_Eta_Docalc_No_Frq_m, only: Spread_Eta_FZP_from_Eta_ZP
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Sparse_Eta_m, only: Sparse_Eta_t

    class(sparse_eta_t), intent(in) :: Eta_FZP_Sparse(:) ! Created by Comp_*_Eta_* above
    real(rp), intent(inout) :: Eta_FZP(:,:) ! Eta_z * Eta_phi for each state
                           ! vector element. Marked inout so that we can clear
                           ! only the nonzeros without making the rest of it
                           ! undefined.
    type(grids_t), intent(in) :: Grids_f
    logical, intent(inout), optional :: Do_Calc_FZP(:,:) ! Do_Calc_FZP(p,v)
                           ! indicates whether there is a contribution for 
                           ! state vector element V at point P on the path.
    integer, intent(inout), optional, target :: NZ_FZP(:,:) ! Nonzeros in eta_FZP
    integer, intent(inout), optional, target :: NNZ_FZP(:)  ! Number of nonzeros in eta_FZP
    integer, intent(in), optional :: One_Sps ! Do extraction for only this species
    class(sparse_eta_t), intent(in), optional :: Eta_ZP_Sparse(:) ! Created by Comp_*_Eta_* above
    logical, intent(in), optional :: Spread  ! For a frequency-dependent species,
                                             ! spread the columns from
                                             ! Eta_ZP_Sparse

    integer :: C, L        ! Column, Last column of Eta_FZP before species S
    logical :: MySpread
    integer :: N_Eta_ZP2   ! Second dimension of Eta_ZP
    integer :: NC          ! Number of columns of Eta_FZP for species S
    integer :: NP          ! Path length, first dimension of Eta_FZP
    integer :: S, S1, S2   ! Species index = index in Eta_FZP_Sparse
    logical :: SpreadIt

    if ( present(one_sps) ) then
      s1 = one_sps
      s2 = one_sps
    else
      s1 = 1
      s2 = size(eta_fzp_sparse,1)
    end if

    mySpread = .false.
    n_eta_zp2 = 0
    if ( present(spread) .and. present(eta_zp_sparse) ) then
      mySpread = spread
      n_eta_zp2 = merge(grids_f%p_len,0,mySpread)
    end if

    np = size(eta_fzp,1)

    block

      real(rp) :: Eta_ZP(np,n_eta_zp2)
      logical :: Do_Calc_ZP(np,n_eta_zp2)
      integer :: NZ_ZP(np,n_eta_zp2)
      integer :: NNZ_ZP(n_eta_zp2)

      do_calc_zp = .false.
      eta_zp = 0
      nnz_zp = 0
      nz_zp = 0
      do s = s1, s2
        if ( .not. allocated(eta_fzp_sparse(s)%rows) .or. &
           & .not. allocated(eta_fzp_sparse(s)%cols) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Rows or Cols of eta_fzp_sparse(s) not allocated in Get_Eta_DoCalc_fzp' )
          return
        end if
        if ( mySpread ) spreadIt =  grids_f%l_f(s-1)+1 /= grids_f%l_f(s)
        if ( spreadIt ) then
          l = grids_f%l_zp(s-1)
          nc = grids_f%l_zp(s)
          do_calc_zp(:,l+1:nc) = .false.
          eta_zp(:,l+1:nc) = 0
          nnz_zp(l+1:nc) = 0
          nz_zp(:,l+1:nc) = 0

          call get_one_eta_docalc_fzp ( eta_zp_sparse(s), eta_zp(:,l+1:nc), &
                                      & do_calc_zp(:,l+1:nc), &
                                      & nz_zp(:,l+1:nc), nnz_zp(l+1:nc) )
          call  spread_eta_fzp_from_eta_zp ( grids_f, eta_zp(:np,:), &
            & do_calc_zp(:np,:), eta_fzp(:np,:), do_calc_fzp(:np,:), &
            & nz_zp(:np,:), nnz_zp, nz_fzp(:np,:), nnz_fzp, s )

        else
          l = grids_f%l_v(s-1)
          nc = grids_f%l_v(s)
          call get_one_eta_docalc_fzp ( eta_fzp_sparse(s), eta_fzp(:,l+1:nc), &
                                      & do_calc_fzp(:,l+1:nc), &
                                      & nz_fzp(:,l+1:nc), nnz_fzp(l+1:nc) )
          if ( present(do_calc_fzp) ) then
            do c = l+1, nc
              do_calc_fzp(:,c) = do_calc_fzp(:,c) .and. &
                               & grids_f%deriv_flags(c)
            end do
          end if
        end if
      end do ! s
 
    end block

  end subroutine Get_Eta_DoCalc_FZP

  subroutine Get_One_Eta_DoCalc_FZP ( Eta_FZP_Sparse, Eta_FZP, &
                                    & Do_Calc_FZP, NZ_FZP, NNZ_FZP )

  ! Compute Eta_FZP and Do_Calc_FZP from Eta_FZP_Sparse for one species.  When the
  ! forward model no longer needs these, this subroutine can be deleted.

    use Get_Eta_Matrix_m, only: Clean_Out_Nonzeros
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    class(sparse_eta_t), intent(in) :: Eta_FZP_Sparse ! Created by Comp_*_Eta_* above
    real(rp), intent(inout) :: Eta_FZP(:,:) ! The columns of Eta_z * Eta_phi
                           ! for one species.  Each column is for one state
                           ! vector element.  Marked inout so that we can clear
                           ! only the nonzeros without making the rest of it
                           ! undefined.
    logical, intent(inout), optional :: Do_Calc_FZP(:,:) ! Do_Calc_FZP(p,v)
                           ! indicates whether there is a contribution for 
                           ! state vector element V at point P on the path.
                           ! Same shape and correspondence for state vector.
    integer, intent(inout), optional, target :: NZ_FZP(:,:) ! Nonzeros in eta_FZP
                           ! Same shape and correspondence for state vector.
    integer, intent(inout), optional, target :: NNZ_FZP(:)  ! Number of nonzeros
                           ! in eta_FZP.  Same columns as Eta_FZP.

    integer :: I, J, K
    integer :: NC 

    nc = size(eta_fzp_sparse%cols)
    call clean_out_nonzeros ( eta_fzp, do_calc_fzp, nz_fzp, nnz_fzp )
    do j = 1, nc
      k = eta_fzp_sparse%cols(j)
      if ( k == 0 ) cycle ! Column is empty
      i = k
      do
        i = eta_fzp_sparse%e(i)%nc ! Next row in this column
        eta_fzp ( eta_fzp_sparse%e(i)%r,j) = eta_fzp_sparse%e(i)%v
        if ( present(nz_fzp) ) then
          nnz_fzp(j) = nnz_fzp(j) + 1
          nz_fzp(nnz_fzp(j),j) = eta_fzp_sparse%e(i)%r
        end if
        if ( present(do_calc_fzp) ) do_calc_fzp ( eta_fzp_sparse%e(i)%r,j ) = .true.
        if ( i == k ) exit        ! No more rows in this column
      end do
    end do

  end subroutine Get_One_Eta_DoCalc_FZP

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
! Revision 2.7  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.6  2018/03/07 17:43:57  pwagner
! Made consistent with new Sparse_t datatype
!
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
