! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Sps_Path_Sparse_m

  implicit NONE

  private
  public :: Comp_Sps_Path_Sparse
  public :: Comp_Sps_Path_Sparse_Frq, Comp_1_Sps_Path_Sparse_Frq
  public :: Comp_Sps_Path_Sparse_No_Frq, Comp_1_Sps_Path_Sparse_No_Frq
  public :: Comp_1_Sps_Path_Sparse_No_Frq_No_Eta
  public :: Comp_1_Sps_Path_Sparse_No_Frq_2D

  interface Comp_Sps_Path_Sparse
    module procedure Comp_Sps_Path_Sparse_Frq, Comp_1_Sps_Path_Sparse_Frq
    module procedure Comp_Sps_Path_Sparse_No_Frq, Comp_1_Sps_Path_Sparse_No_Frq
    module procedure Comp_1_Sps_Path_Sparse_No_Frq_No_Eta
    module procedure Comp_1_Sps_Path_Sparse_No_Frq_2D
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Comp_Sps_Path_Sparse_Frq ( Grids_f, Frq, Eta_ZP, Eta_FZP, &
                                      & Sps_Path, LO, Sideband )

    ! Compute the Sps_Path for species that are frequency dependent.
    ! This assumes that it has already been computed for species that
    ! are not frequency dependent.

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP, R8
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f        ! Quantity values
    real(r8), intent(in) :: Frq                 ! Frequency at which to compute
                                                ! values in Sps_Path.
    type(sparse_eta_t), intent(in) :: Eta_ZP(:) ! Interpolate Zeta X H
                                                ! to Sps_Path, same size as
                                                ! Grids_f%Mol
    type(sparse_eta_t), intent(inout) :: Eta_FZP(:) ! Interpolate
                                                ! F X Zeta X H to Sps_Path, size
                                                ! is number of quantities in
                                                ! Grids_f that have more than one
                                                ! element in Frq_Basis.
    real(rp), intent(inout) :: Sps_Path(:,:)    ! Path X Sps -- VMR values.
    real(r8), intent(in) :: LO                  ! Local oscillator frequency, GHz
    integer, intent(in) :: Sideband             ! -1, 1, or 0.  Zero means
                                                ! quantities' frequency bases
                                                ! absolute, not I.F.

    integer :: I

    do i = 1, size(eta_fzp)
      if ( grids_f%l_f(i-1)+1 == grids_f%l_f(i) ) cycle ! No frequency dependency for this species
      call comp_1_sps_path_sparse_frq ( grids_f, Frq, Eta_ZP(i), Eta_FZP(i), &
                                      & Sps_Path(:,i), LO, Sideband, i )
    end do

  end subroutine Comp_Sps_Path_Sparse_Frq

  subroutine Comp_1_Sps_Path_Sparse_Frq ( Grids_f, Frq, Eta_ZP, Eta_FZP, &
                                        & Sps_Path, LO, Sideband, N )

    ! Compute the Sps_Path for species that are frequency dependent.
    ! This assumes that it has already been computed for species that
    ! are not frequency dependent.

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP, R8
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: grids_f     ! Quantity values
    real(r8), intent(in) :: Frq              ! Frequency at which to compute
                                             ! values in Sps_Path.
    type(sparse_eta_t), intent(in) :: Eta_ZP ! Interpolate Zeta X H to Sps_Path.
    type(sparse_eta_t), intent(inout) :: Eta_FZP ! Interpolate F X Zeta X H to
                                             ! Sps_Path.
    real(rp), intent(inout) :: Sps_Path(:)   ! Path for 1 Sps -- VMR values.
    real(r8), intent(in) :: LO               ! Local oscillator frequency, GHz
    integer, intent(in) :: Sideband          ! -1, 1, or 0.  Zero means
                                             ! quantities' frequency bases
                                             ! absolute, not I.F.
    integer, intent(in), optional :: N       ! Which quantity, default 1

    type(sparse_eta_t) :: Eta_F
    integer :: F1, F2
    integer :: MyN
    integer :: What

    myN = 1
    if ( present(n) ) myN = n

    what = grids_f%qtyStuff(myN)%qty%template%name

    ! Compute Eta_F for quantity N.  We don't need it for anything other than
    ! computing Eta_FZP.
    f1 = grids_f%l_f(myN-1)+1
    f2 = grids_f%l_f(myN)
    ! Eta_1D isn't prepared to work with a basis that's not in increasing
    ! order.  lo - grids_f%frq_basis(f1:f2) would be in decreasing order.
    ! lo - grids_f%frq_basis(f2:f1:-1) would be in increasing order, but the
    ! column subscripts in Eta would be inverted.  So for the lower sideband,
    ! we use -(lo-grids_f%frq_basis(f1:f2)) and -Frq, which produces the
    ! correct coefficients, in the correct order.

    ! The brackets around the first argument aren't strictly necessary because
    ! the first term is an array, but without them ifort 17 appears to believe
    ! the result is a zero-size array, so Eta_f%Eta_1D doesn't do anything.
!! Even with brackets, ifort 17 doesn't work, but this apparent opportunity
!! to use grids_f%frq_basis(f1:f2)+sideband*lo as the actual argument to another
!! subroutine convinces it to do the right thing.  The call doesn't ever happen
!! because array dimension extents are never negative.  If ifort 18 works....
if ( size(sps_path,1) < 0 ) print '(1p5g15.6)', grids_f%frq_basis(f1:f2)+sideband*lo
    call eta_f%eta_1d ( [ grids_f%frq_basis(f1:f2) + sideband*lo ], &
                      & [ merge(frq, sideband*frq, sideband==0) ], what=what )
!! and this one appears to be necessary too
if ( size(sps_path,1) < 0 ) call eta_f%dump ( name='Eta_F', width=4 )
    call eta_fzp%eta_nd ( eta_f, eta_zp, what=what, resize=.true., one_row_ok=.true. )

    ! Now that we have Eta_FZP, we can finally interpolate.
    ! Sps_path = Eta_fzp .dot. Grids_f%c(myN)%v1
    call eta_fzp%sparse_dot_vec ( grids_f%c(myN)%v1, sps_path )
    if ( grids_f%lin_log(myN) ) sps_path = exp(sps_path)

  end subroutine Comp_1_Sps_Path_Sparse_Frq

  subroutine Comp_Sps_Path_Sparse_No_Frq ( Grids_f, Eta_ZP, Sps_Path, Eta_FZP )

    ! Compute the Sps_Path for species that are not frequency dependent.
    ! Compute the indices between Eta_ZP and Eta_FZP for frequency-dependent
    ! quantities.

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f           ! Quantity values
    type(sparse_eta_t), intent(inout) :: Eta_ZP(:) ! Interpolate Zeta X H
                                                   ! to Sps_Path, same size as
                                                   ! Grids_f%Mol.
    real(rp), intent(inout) :: Sps_Path(:,:)       ! Path X Sps -- VMR values.
    type(sparse_eta_t), intent(inout), optional :: Eta_FZP(:)  ! Interpolate F X
                                                   ! Zeta X H to Sps_Path, size
                                                   ! is the number of quantities
                                                   ! in Grids_f that have more
                                                   ! than one element in
                                                   ! Frq_Basis.

    integer :: Sps

    if ( present(eta_fzp) ) then
      do sps = 1, size(eta_zp)
        if ( grids_f%l_f(sps-1)+1 == grids_f%l_f(sps) ) & ! Not frequency dependent
          ! Sps_path(:,sps) = Eta_zp(sps) .dot. Grids_f%c(sps)%v1
          & call eta_fzp(sps)%sparse_dot_vec ( grids_f%c(sps)%v1, sps_path(:,sps) )
        if ( grids_f%lin_log(sps) ) sps_path(:,sps) = exp(sps_path(:,sps))
      end do
    else
      do sps = 1, size(eta_zp)
        call eta_zp(sps)%sparse_dot_vec ( grids_f%c(sps)%v1, sps_path(:,sps) )
        if ( grids_f%lin_log(sps) ) sps_path(:,sps) = exp(sps_path(:,sps))
      end do
    end if

  end subroutine Comp_Sps_Path_Sparse_No_Frq

  subroutine Comp_1_Sps_Path_Sparse_No_Frq ( Grids_f, Eta_ZP, Sps_Path, N )

    ! Compute the Sps_Path for one species that is not frequency dependent.

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f    ! Quantity values
    type(sparse_eta_t), intent(inout) :: Eta_ZP ! Interpolate Zeta X H to Sps_Path
    real(rp), intent(inout) :: Sps_Path(:)  ! Path for 1 Sps -- VMR values.
    integer, intent(in), optional :: N      ! Which species, default 1

    integer :: MyN

    myN = 1
    if ( present(n) ) myN = n

    call eta_zp%sparse_dot_vec ( grids_f%c(myN)%v1, sps_path )
    if ( grids_f%lin_log(myN) ) sps_path = exp(sps_path)

  end subroutine Comp_1_Sps_Path_Sparse_No_Frq

  subroutine Comp_1_Sps_Path_Sparse_No_Frq_No_Eta ( Grids_f, Tan_Pt, Z_Path, &
                                                  & Phi_Path, Sps_Path, N )

    ! Compute the Sps_Path for one species that is not frequency dependent.

    use Comp_Eta_Docalc_Sparse_m, only: Comp_Eta
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f    ! Quantity values
    integer, intent(in) :: Tan_Pt           ! To split Z_Path into two
                                            ! monotone halves
    real(rp), intent(in) :: Z_Path(:)       ! Zetas on the path
    real(rp), intent(in) :: Phi_Path(:)     ! Phis on the path, Radians
    real(rp), intent(inout) :: Sps_Path(:)  ! Path for 1 Sps -- VMR values.
    integer, intent(in), optional :: N      ! Which species, default 1

    type(sparse_eta_t) :: Eta_ZP            ! Interpolate Zeta X H to Sps_Path
    integer :: MyN

    myN = 1
    if ( present(n) ) myN = n

    call comp_eta ( grids_f, tan_pt, z_path, phi_path, eta_zp )
    call eta_zp%sparse_dot_vec ( grids_f%c(myN)%v1, sps_path )
    if ( grids_f%lin_log(myN) ) sps_path = exp(sps_path)

  end subroutine Comp_1_Sps_Path_Sparse_No_Frq_No_Eta

  subroutine Comp_1_Sps_Path_Sparse_No_Frq_2D ( Grids_f, Eta_ZP, Sps_Path, N )

    ! Compute the Sps_Path for all the columns of one species that is not
    ! frequency dependent, but is represented as if it were, e.g. magnetic
    ! field, for which Sps_Path has three columns, but Grids_f%c(n)%v4 has
    ! three ROWS.

    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    type(grids_t), intent(in) :: Grids_f     ! Quantity values
    type(sparse_eta_t), intent(inout) :: Eta_ZP ! Interpolate Zeta X H to Sps_Path,
    real(rp), intent(inout) :: Sps_Path(:,:) ! Path for 1 Sps, e.g. magnetic field
    integer, intent(in), optional :: N       ! Which species, default 1

    integer :: I, MyN, NC

    myN = 1
    if ( present(n) ) myN = n

    nc = grids_f%l_f(myN) - grids_f%l_f(myN-1)
    ! v4 is Frequency X Zeta X Phi.  We want to compute the product
    ! Eta * v4(i,:,:,1) and put it in sps_path(:,i).  We can't take a 1-D
    ! pointer to v4(i,:,:) because it's not contiguous.
    do i = 1, size(sps_path,2) ! should be the same as size(grids_f%c(myN)%v4,1)
      call eta_zp%sparse_dot_vec_2d ( grids_f%c(myN)%v4(i,:,:,1), sps_path(:,i) )
    end do
    if ( grids_f%lin_log(myN) ) sps_path = exp(sps_path)

  end subroutine Comp_1_Sps_Path_Sparse_No_Frq_2D

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Comp_Sps_Path_Sparse_m

! $Log$
! Revision 2.7  2018/09/12 23:50:36  vsnyder
! Make N argument of Comp_1_Sps_Path_Sparse_Frq optional
!
! Revision 2.6  2018/09/12 22:02:13  vsnyder
! Use Comp_Eta instead of Comp_Eta_DoCalc_Sparse
!
! Revision 2.5  2018/09/05 20:54:58  vsnyder
! Add Comp_1_Sps_Path_Sparse_No_Frq_No_Eta and Comp_1_Sps_Path_Sparse_No_Frq_2D
!
! Revision 2.4  2018/08/28 22:15:04  vsnyder
! Make Eta_FZP optional in Comp_Sps_Path_Sparse_No_Frq
!
! Revision 2.3  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.2  2017/03/11 00:53:05  vsnyder
! Use Grids_f instead of Qty_Stuff, remove Lists_F_t types
!
! Revision 2.1  2017/01/17 19:57:18  vsnyder
! Initial commit
!
