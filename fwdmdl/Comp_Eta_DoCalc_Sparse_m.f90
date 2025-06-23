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

  public :: Comp_Eta, Comp_Eta_Docalc_Sparse
  public :: Comp_All_Eta_2D,  Comp_All_Eta_2D_Local_1D, Comp_One_Eta_2D
  public :: Comp_All_Eta_QTM, Comp_One_Eta_QTM
  public :: Comp_One_Eta_Z

  interface Comp_Eta
    module procedure Comp_All_Eta_2D, Comp_All_Eta_2D_Local_1D, Comp_One_Eta_2D
    module procedure Comp_One_Eta_2D_Local_1D
    module procedure Comp_All_Eta_QTM, Comp_One_Eta_QTM
  end interface

  interface Comp_Eta_Docalc_Sparse
    module procedure Comp_All_Eta_2D, Comp_All_Eta_2D_Local_1D, Comp_One_Eta_2D
    module procedure Comp_One_Eta_2D_Local_1D
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

  ! Compute interpolation coefficients for Zeta x Phi for all molecules,
  ! assuming some will be frequency dependent, and it is necessary to
  ! keep Eta_Phi and Eta_Z.  If those assumptions are not necessary, use
  ! Comp_All_Eta_2D_Local_1D.

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
    class(sparse_eta_t), intent(inout), optional :: Eta_fzp(:) ! Created here
    integer, intent(in), optional :: Skip      ! At tangent point, default 1

    integer :: I

    if ( present(eta_fzp) ) then
      do i = 1, size(grids_f%mol)
        if ( grids_f%l_f(i-1)+1 == grids_f%l_f(i) ) then ! Not frequency dependent
          call comp_one_eta_2d ( grids_f, tan_pt, z_path, eta_z(i), &
                               & phi_path, eta_phi(i), eta_fzp(i), i, skip )
        else ! Frequency dependent, fill Eta_ZP so we can calculate Eta_FZP
             ! when we have frequency
          call comp_one_eta_2d ( grids_f, tan_pt, z_path, eta_z(i), &
                               & phi_path, eta_phi(i), eta_zp(i), i, skip )
        end if
      end do
    else
      call comp_one_eta_2d ( grids_f, tan_pt, z_path, eta_z(i), &
                           & phi_path, eta_phi(i), eta_zp(i), i, skip )
    end if

  end subroutine Comp_All_Eta_2D

  ! -----------------------------------  Comp_All_Eta_2D_Local_1D  -----

  subroutine Comp_All_Eta_2D_Local_1D ( Grids_f, Tan_Pt, Z_Path, Phi_Path, &
                                      & Eta_zp, Skip )

  ! Compute interpolation coefficients for Zeta x Phi for all molecules
  ! without provision for frequency dependence.  Use this one if it isn't
  ! necessary to keep Eta_Phi and Eta_Z

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f       ! Quantity values
    integer, intent(in) :: Tan_Pt              ! To split Z_Path into two
                                               ! monotone halves
    real(rp), intent(in) :: Z_Path(:)          ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)        ! Phis on the path
    class(sparse_eta_t), intent(inout) :: Eta_zp(:)  ! Created here
    integer, intent(in), optional :: Skip      ! At tangent point, default 1

    integer :: I

    do i = 1, size(grids_f%l_v)
      call comp_one_eta_2d_local_1D ( grids_f, tan_pt, z_path, phi_path, &
                                    & eta_zp(i), i, skip )
    end do

  end subroutine Comp_All_Eta_2D_Local_1D

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
      call comp_one_eta_QTM ( grids_f, tan_pt, z_path, eta_z(i), &
                            & eta_path(i), eta_zQ(i), i, skip )
    end do

  end subroutine Comp_All_Eta_QTM

  ! --------------------------------------------  Comp_One_Eta_2D  -----

  subroutine Comp_One_Eta_2D ( Grids_f, Tan_Pt, Z_Path, Eta_Z, &
                             & Phi_Path, Eta_Phi, Eta_zP, N, Skip )

  ! Compute interpolation coefficients for Zeta x Phi for one molecule

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f         ! Quantity values
    integer, intent(in) :: Tan_Pt                ! To split Z_Path into two
                                                 ! monotone halves
    real(rp), intent(in) :: Z_Path(:)            ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z    ! from VMR's Zeta to path.
    real(rp), intent(in) :: Phi_Path(:)          ! Phis on the path, Radians
    class(sparse_eta_t), intent(out) :: Eta_Phi  ! from VMR's Phi to path.
    class(sparse_eta_t), intent(inout) :: Eta_zP ! Created here
    integer, intent(in), optional :: N           ! Which quantity, default 1
    integer, intent(in), optional :: Skip        ! At tangent point, default 1

    integer :: MyN
    integer :: P1, P2 ! Boundaries from Grids_f%l_p
    integer :: What

    myN = 1
    if ( present(n) ) myN = n

    what = grids_f%qtyStuff(myN)%qty%template%name

    p1 = grids_f%l_p(myN-1)+1
    p2 = grids_f%l_p(myN)

    ! Create and compute Eta_Phi
    call eta_phi%eta_1d ( grids_f%phi_basis(p1:p2), phi_path, &
                        & what=what, resize=.true. )
    eta_phi%lbnd = [ p1 ] ! Column bounds
    eta_phi%ubnd = [ p2 ] !   in phi_basis

    ! Create Eta_Z
    call comp_one_eta_z ( grids_f, tan_pt, z_path, eta_z, myN, skip )

    ! Compute Eta_ZP
    call eta_zp%eta_nd ( eta_z, eta_phi, what=what, resize=.true. )

  end subroutine Comp_One_Eta_2D

  ! -----------------------------------  Comp_One_Eta_2D_Local_1D  -----

  subroutine Comp_One_Eta_2D_Local_1D ( Grids_f, Tan_Pt, Z_Path, &
                                      & Phi_Path, Eta_zP, N, Skip )

  ! Compute interpolation coefficients for Zeta x Phi for one molecule.
  ! Use this one if you don't need to keep Eta_Z and Eta_Phi

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f         ! Quantity values
    integer, intent(in) :: Tan_Pt                ! To split Z_Path into two
                                                 ! monotone halves
    real(rp), intent(in) :: Z_Path(:)            ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)          ! Phis on the path, Radians
    class(sparse_eta_t), intent(inout) :: Eta_zP ! Created here
    integer, intent(in), optional :: N           ! Which quantity, default 1
    integer, intent(in), optional :: Skip        ! At tangent point, default 1

    type(sparse_eta_t) :: Eta_Z    ! from VMR's Zeta to path.
    type(sparse_eta_t) :: Eta_Phi  ! from VMR's Phi to path.
    integer :: MyN
    integer :: P1, P2 ! Boundaries from Grids_f%l_p
    integer :: What

    myN = 1
    if ( present(n) ) myN = n

    what = grids_f%qtyStuff(myN)%qty%template%name

    p1 = grids_f%l_p(myN-1)+1
    p2 = grids_f%l_p(myN)

    ! Create and compute Eta_Phi
    call eta_phi%eta_1d ( grids_f%phi_basis(p1:p2), phi_path, &
                        & what=what, resize=.true. )
    eta_phi%lbnd = [ p1 ] ! Column bounds
    eta_phi%ubnd = [ p2 ] !   in phi_basis

    ! Create Eta_Z
    call comp_one_eta_z ( grids_f, tan_pt, z_path, eta_z, myN, skip )

    ! Compute Eta_ZP
    call eta_zp%eta_nd ( eta_z, eta_phi, what=what, resize=.true. )

  end subroutine Comp_One_Eta_2D_Local_1D

  ! -------------------------------------------  Comp_One_Eta_QTM  -----

  subroutine Comp_One_Eta_QTM ( Grids_f, Tan_Pt, Z_Path, Eta_Z, Eta_Path, &
                              & Eta_zQ, N, Skip )

  ! Compute interpolation coefficients for Zeta x QTM for one molecule

    use Sparse_Eta_m, only: Sparse_Eta_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP

    type(grids_t), intent(in) :: Grids_f        ! Quantity values
    integer, intent(in) :: Tan_Pt               ! To split Z_Path into two
                                                ! monotone halves
    real(rp), intent(in) :: Z_Path(:)           ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z   ! from VMR's zeta to path.
    class(sparse_eta_t), intent(in) :: Eta_Path ! from QTM to path.  Not
                                                ! computed here because there's
                                                ! only one QTM for everything
    class(sparse_eta_t), intent(out) :: Eta_zQ  ! Created here
    integer, intent(in), optional :: N          ! Which quantity, default 1
    integer, intent(in), optional :: Skip       ! At tangent point, default 1

    ! Create Eta_Z
    call comp_one_eta_z ( grids_f, tan_pt, z_path, eta_z, n, skip )

    ! Compute Eta_ZQ
!     call eta_zp%eta_nd ( eta_z, eta_path, what=what, resize=.true. )

  end subroutine Comp_One_Eta_QTM

  ! ---------------------------------------------  Comp_One_Eta_Z  -----

  subroutine Comp_One_Eta_Z ( Grids_f, Tan_Pt, Z_Path, Eta_Z, &
                            & N, Skip )

  ! Compute interpolation coefficients for Zeta for one molecule

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f         ! Quantity values
    integer, intent(in) :: Tan_Pt                ! To split Z_Path into two
                                                 ! monotone halves
    real(rp), intent(in) :: Z_Path(:)            ! Zetas on the path
    class(sparse_eta_t), intent(out) :: Eta_Z    ! from VMR's Zeta to path.
    integer, intent(in), optional :: N           ! Which quantity, default 1
    integer, intent(in), optional :: Skip        ! At tangent point, default 1

    integer :: Z1, Z2 ! Boundaries from Grids_f%l_z
    integer :: MyN, MySkip
    integer :: N_Path
    integer :: What

    myN = 1
    if ( present(n) ) myN = n

    mySkip = 1
    if ( present(skip) ) mySkip = skip

    what = grids_f%qtyStuff(myN)%qty%template%name

    n_path = size(z_path) ! Assumed == size(p_path)

    z1 = grids_f%l_z(myN-1)+1
    z2 = grids_f%l_z(myN)

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
! Revision 2.10  2018/09/12 22:01:15  vsnyder
! Add Comp_Eta generic, identical to Comp_Eta_DoCalc_Sparse.  Add
! Comp_All_Eta_2D_Local_1D.
!
! Revision 2.9  2018/08/28 22:14:02  vsnyder
! Add Comp_One_Eta_2D_Local_1D, make some arguments optional
!
! Revision 2.8  2018/05/17 01:31:52  vsnyder
! Remove no-longer-used routines
!
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
