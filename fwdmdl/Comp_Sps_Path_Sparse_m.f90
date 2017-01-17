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

  use Indexed_Values_m, only: Value_2D_Lists_t, Value_3D_Lists_t, &
    & Value_QTM_2D_Lists_t

  implicit NONE

  private
  public :: Comp_Sps_Path_Sparse
  public :: Comp_Sps_Path_Sparse_Frq, Comp_1_Sps_Path_Sparse_Frq
  public :: Comp_Sps_Path_Sparse_No_Frq, Comp_1_Sps_Path_Sparse_No_Frq
  public :: Get_Frequency_Index_For_Eta_FZP

  interface Comp_Sps_Path_Sparse
    module procedure Comp_Sps_Path_Sparse_Frq, Comp_1_Sps_Path_Sparse_Frq
    module procedure Comp_Sps_Path_Sparse_No_Frq, Comp_1_Sps_Path_Sparse_No_Frq
  end interface

  type, public, extends(value_2d_lists_t) :: Value_2d_Lists_F_t
  ! integer :: N                             ! Inherited
  ! type(value_2d_lists_t), allocatable :: Eta(:) ! Inherited
    integer :: Frq_Index = 0 ! Index in the "other" Eta_[F]ZP
  end type Value_2d_Lists_F_t

  type, public, extends(value_3d_lists_t) :: Value_3d_Lists_F_t
  ! integer :: N                             ! Inherited
  ! type(value_2d_lists_t), allocatable :: Eta(:) ! Inherited
    integer :: Frq_Index = 0 ! Index in the "other" Eta_[F]ZP
  end type Value_3d_Lists_F_t

  type, public, extends(value_QTM_2d_lists_t) :: Value_QTM_2D_Lists_f_t
  ! integer :: N                             ! Inherited
  ! type(value_QTM_2d_lists_t), allocatable :: Eta(:) ! Inherited
    integer :: Frq_Index = 0 ! Index in the "other" Eta_[F]ZP
  end type Value_QTM_2D_Lists_f_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Comp_Sps_Path_Sparse_Frq ( Qty_Stuff, Frq, Eta_ZP, Eta_FZP, &
                                      & Sps_Path, LO, Sideband )

    ! Compute the Sps_Path for species that are frequency dependent.
    ! This assumes that it has already been computed for species that
    ! are not frequency dependent.

    use ForwardModelConfig, only: QtyStuff_t
    use MLSKinds, only: RP, R8

    type(qtyStuff_t), intent(in) :: Qty_Stuff(:) ! Quantity values on profiles
    real(r8), intent(in) :: Frq                  ! Frequency at which to compute
                                                 ! values in Sps_Path.
    type(value_2d_lists_f_t), intent(in) :: Eta_ZP(:) ! Interpolate Zeta X H
                                                 ! to Sps_Path, same size as
                                                 ! Qty_Stuff.
    type(value_3d_lists_f_t), intent(inout) :: Eta_FZP(:) ! Interpolate
                                                 ! F X Zeta X H to Sps_Path, size
                                                 ! is number of elements of
                                                 ! Qty_Stuff that have a
                                                 ! Frequencies component.
    real(rp), intent(inout) :: Sps_Path(:,:)     ! Path X Sps -- VMR values.
    real(r8), intent(in) :: LO                   ! Local oscillator frequency, GHz
    integer, intent(in) :: Sideband              ! -1, 1, or 0.  Zero means
                                                 ! quantities' frequency bases
                                                 ! absolute, not I.F.

    integer :: I, Q

    do i = 1, size(eta_fzp)
      q = eta_fzp(i)%frq_index ! index in Eta_ZP, Qty_Stuff, and Sps_Path(:,q)
      call comp_1_sps_path_sparse_frq ( Qty_Stuff(q), Frq, Eta_ZP, Eta_FZP(i), &
                                      & Sps_Path(:,q), LO, Sideband )
    end do

  end subroutine Comp_Sps_Path_Sparse_Frq

  subroutine Comp_1_Sps_Path_Sparse_Frq ( Qty_Stuff, Frq, Eta_ZP, Eta_FZP, &
                                        & Sps_Path, LO, Sideband )

    ! Compute the Sps_Path for species that are frequency dependent.
    ! This assumes that it has already been computed for species that
    ! are not frequency dependent.

    use ForwardModelConfig, only: QtyStuff_t
    use Get_Eta_List_m, only: Get_Eta_List
    use Indexed_Values_m, only: Interpolate, Value_1D_List_t
    use MLSKinds, only: RP, R8

    type(qtyStuff_t), intent(in) :: Qty_Stuff    ! Quantity values on profiles
    real(r8), intent(in) :: Frq                  ! Frequency at which to compute
                                                 ! values in Sps_Path.
    type(value_2d_lists_f_t), intent(in) :: Eta_ZP(:) ! Interpolate Zeta X H
                                                 ! to Sps_Path.
    type(value_3d_lists_f_t), intent(inout) :: Eta_FZP ! Interpolate F X Zeta X H
                                                 ! to Sps_Path.
    real(rp), intent(inout) :: Sps_Path(:)       ! Path for 1 Sps -- VMR values.
    real(r8), intent(in) :: LO                   ! Local oscillator frequency, GHz
    integer, intent(in) :: Sideband              ! -1, 1, or 0.  Zero means
                                                 ! quantities' frequency bases
                                                 ! absolute, not I.F.

    type(value_1D_list_t) :: Eta_F(1)
    integer :: N_F, Q

    q = eta_fzp%frq_index
    ! Compute Eta_F for quantity.  We don't need it for anything else.
    n_f = size(qty_stuff%qty%template%frequencies)
    select case ( sideband )
    case ( -1 )
      call get_eta_list ( lo-qty_stuff%qty%template%frequencies(n_f:1:-1), &
                        & [ frq ], eta_f )
    case ( +1 )
      call get_eta_list ( lo+qty_stuff%qty%template%frequencies(1:n_f:+1), &
                        & [ frq ], eta_f )
    case ( 0 )
      call get_eta_list ( qty_stuff%qty%template%frequencies, &
                        & [ frq ], eta_f )
    end select

    call get_eta_list ( eta_f, eta_zp(q)%eta(:eta_zp(q)%n), &
                      & eta_fzp%eta(:eta_fzp%n) )

    ! Now that we have Eta_FZP, we can finally interpolate.
    call interpolate ( qty_stuff%qty%value3, eta_fzp%value_3d_lists_t, sps_path )

  end subroutine Comp_1_Sps_Path_Sparse_Frq

  subroutine Comp_Sps_Path_Sparse_No_Frq ( Qty_Stuff, Eta_ZP, Eta_FZP, Sps_Path )

    ! Compute the Sps_Path for species that are not frequency dependent.
    ! Compute the indices between Eta_ZP and Eta_FZP for frequency-dependent
    ! quantities.

    use ForwardModelConfig, only: QtyStuff_t
    use Indexed_Values_m, only: Interpolate
    use MLSKinds, only: RP

    type(qtyStuff_t), intent(in) :: Qty_Stuff(:) ! Quantity values on profiles
    type(value_2d_lists_f_t), intent(inout) :: Eta_ZP(:) ! Interpolate Zeta X H
                                                 ! to Sps_Path, same size as
                                                 ! Qty_Stuff.
    type(value_3d_lists_f_t), intent(inout) :: Eta_FZP(:) ! Interpolate F X Zeta X H
                                                 ! to Sps_Path, size is number
                                                 ! of elements of Qty_Stuff that
                                                 ! have a Frequencies component.
    real(rp), intent(inout) :: Sps_Path(:,:)     ! Path X Sps -- VMR values.

    integer :: I

    call get_frequency_index_for_eta_fzp ( qty_stuff, eta_zp, eta_fzp )
    do i = 1, size(eta_zp)
      if ( eta_zp(i)%frq_index == 0 ) &
        & call interpolate ( qty_stuff(i)%qty%values, eta_zp(i)%value_2d_lists_t, &
                           & sps_path(:,i) )
    end do

  end subroutine Comp_Sps_Path_Sparse_No_Frq

  subroutine Comp_1_Sps_Path_Sparse_No_Frq ( Qty_Stuff, Eta_ZP, Sps_Path )

    ! Compute the Sps_Path for species that are not frequency dependent.
    ! Compute the indices between Eta_ZP and Eta_FZP for frequency-dependent
    ! quantities.

    use ForwardModelConfig, only: QtyStuff_t
    use Indexed_Values_m, only: Interpolate => Interpolate_Polymorphic
    use MLSKinds, only: RP

    type(qtyStuff_t), intent(in) :: Qty_Stuff    ! Quantity values on profiles
    type(value_2d_lists_f_t), intent(inout) :: Eta_ZP ! Interpolate Zeta X H
                                                 ! to Sps_Path, same size as
                                                 ! Qty_Stuff.
    real(rp), intent(inout) :: Sps_Path(:)       ! Path for 1 Sps -- VMR values.

    if ( eta_zp%frq_index == 0 ) &
      & call interpolate ( qty_stuff%qty%values, eta_zp%eta, sps_path )

  end subroutine Comp_1_Sps_Path_Sparse_No_Frq

  subroutine Get_Frequency_Index_For_Eta_FZP ( Qty_Stuff, Eta_ZP, Eta_FZP )
    ! Compute the index in Eta_FZP for elements of Qty_Stuff that
    ! have an allocated Frequencies component of nonzero size.
    use ForwardModelConfig, only: QtyStuff_T
    type(qtyStuff_t), intent(in) :: Qty_Stuff(:) ! Quantity values on profiles
    type(value_2d_lists_f_t), intent(inout) :: Eta_ZP(:)   ! Interpolator
    type(value_3d_lists_f_t), intent(inout) :: Eta_FZP(:)  ! Interpolate
                                                 ! F X Zeta X H to Sps_Path,
                                                 ! size is number of elements of
                                                 ! Qty_Stuff that have a
                                                 ! Frequencies component.

    integer :: FX          ! Frequency index, to put into Eta_ZP(i) to access
                           ! the appropriate element of Eta_FZP
    integer :: I

    fx = 0
    do i = 1, size(qty_stuff)
      eta_zp(i)%frq_index = 0
      if ( associated(qty_stuff(i)%qty%template%frequencies) ) then
        if ( size(qty_stuff(i)%qty%template%frequencies) /= 0 ) then
          fx = fx + i
          eta_zp(i)%frq_index = fx
          eta_fzp(fx)%frq_index = i
        end if
      end if
    end do
  end subroutine Get_Frequency_Index_For_Eta_FZP

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
! Revision 2.1  2017/01/17 19:57:18  vsnyder
! Initial commit
!
