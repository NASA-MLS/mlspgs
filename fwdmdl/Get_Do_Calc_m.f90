! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Do_Calc_m

  implicit NONE
  private

  public :: Clean_Out_Nonzeros, Get_Do_Calc, Get_Eta_Do_Calc

  ! Set nonzeroes into a logical array at places given by subscripts in
  ! an interpolation coefficient list.  This is hopefully a temporary
  ! code until the procedures that care about Do_Calc variables use the
  ! subscripts in interpolation coefficient lists directly.

  interface Get_Do_Calc
    module procedure Get_Do_Calc_2
  end interface

  ! Clear old nonzeros and true values, then put in new ones from a list.
  ! This is hopefully a temporary code until the procedures that care
  ! about Do_Calc variables use the subscripts in interpolation
  ! coefficient lists directly.

  interface Get_Eta_Do_Calc
    module procedure Get_Eta_Do_Calc_2, Get_Eta_Do_Calc_3
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Clean_Out_Nonzeros ( Eta_Array, Do_Calc, NZ, NNZ )
    ! Clean out old nonzeros in Eta_Array and old true values in Do_Calc
    use MLSKinds, only: RP
    real(rp), intent(inout) :: Eta_Array(:,:)
    logical, intent(inout) :: Do_Calc(:,:)      ! Same shape as Eta_Array
    integer, intent(inout) :: NZ(:,:) ! Row subscripts of nonzeros
    integer, intent(inout) :: NNZ(:)  ! How many nonzeros in a column
    integer :: I
    do i = 1, size(eta_array,2)
      eta_array(nz(:nnz(i),i),i) = 0
      do_calc(nz(:nnz(i),i),i) = .false.
      nnz(i) = 0
    end do
  end subroutine Clean_Out_Nonzeros

  subroutine Get_Do_Calc_2 ( Eta, Two_D_Bounds, Do_Calc, NZ, NNZ )
    use Array_Stuff, only: Element_Position
    use Indexed_Values_m, only: Value_2D_List_t

    type(value_2d_list_t), intent(in) :: Eta(:) ! Interpolation coefficient list
    integer, intent(in) :: Two_D_Bounds(2)      ! Zeta and Phi bounds
    logical, intent(out) :: Do_Calc(:,:)        ! Actually 3D; first dimension
      ! is path length, same as Eta(:).  Second dimension is Zeta X Phi for
      ! the state vector.
    integer, intent(inout), optional :: NZ(:,:) ! Row subscripts of nonzeros
    integer, intent(inout), optional :: NNZ(:)  ! How many nonzeros in a column

    integer :: I, J, K

    if ( present(nz) ) then
      ! Clean out old true values in Do_Calc
      do i = 1, size(do_calc,2)
        do_calc(nz(:nnz(i),i),i) = .false.
        nnz(i) = 0
      end do
      do i = 1, size(eta,1)
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], two_d_bounds)
          do_calc(i,k) = .true.
          nnz(k) = nnz(k) + 1
          nz(nnz(k),k) = i
        end do
      end do
    else
      do_calc = .false.
      do i = 1, size(eta,1)
        do j = 1, eta(i)%n
          do_calc( i,element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], &
                                     & two_d_bounds) ) = .true.
        end do
      end do
    end if

  end subroutine Get_Do_Calc_2

  subroutine Get_Eta_Do_Calc_2 ( Eta, Two_D_Bounds, Eta_Array, Do_Calc, NZ, NNZ, &
                               & DerivFlags )
    use Array_Stuff, only: Element_Position
    use Indexed_Values_m, only: Value_2D_List_t
    use MLSKinds, only: RP
    type(value_2d_list_t), intent(in) :: Eta(:) ! Interpolation coefficient list
    integer, intent(in) :: Two_D_Bounds(2)      ! Zeta and Phi bounds
    real(rp), intent(inout) :: Eta_Array(:,:)   ! Actually 3D; first dimension
      ! is path length, same as Eta(:).  Second dimension is Zeta X Phi for
      ! the state vector.
    logical, intent(inout) :: Do_Calc(:,:)      ! Same shape as Eta_Array
    integer, intent(inout), optional :: NZ(:,:) ! Row subscripts of nonzeros
    integer, intent(inout), optional :: NNZ(:)  ! How many nonzeros in a column
    logical, intent(in), optional :: DerivFlags(:) ! Compute derivatives w.r.t.
                                                ! these state-vector elements.

    integer :: I, J, K
    integer :: N_Path

    n_path = size(eta)

    if ( present(nz) ) then
      call clean_out_nonzeros ( eta_array, do_calc, nz, nnz )
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], two_d_bounds)
          eta_array(i,k) = eta(i)%v(j)%v
          do_calc(i,k) = .true.
          nnz(k) = nnz(k) + 1
          nz(nnz(k),k) = i
        end do
      end do
    else
      eta_array = 0
      do_calc = .false.
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], two_d_bounds)
          eta_array(i,k) = eta(i)%v(j)%v
          do_calc(i,k) = .true.
        end do
      end do
    end if

    if ( present(derivFlags) ) then
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], two_d_bounds)
          do_calc(i,k) = do_calc(i,k) .and. derivFlags(k)
        end do
      end do
    end if

  end subroutine Get_Eta_Do_Calc_2

  subroutine Get_Eta_Do_Calc_3 ( Eta, Three_D_Bounds, Eta_Array, Do_Calc, NZ, NNZ, &
                               & DerivFlags )
    use Array_Stuff, only: Element_Position
    use Indexed_Values_m, only: Value_3D_List_t
    use MLSKinds, only: RP
    type(value_3d_list_t), intent(in) :: Eta(:) ! Interpolation coefficient list
    integer, intent(in) :: Three_D_Bounds(3)    ! Freq, Zeta and Phi bounds
    real(rp), intent(inout) :: Eta_Array(:,:)   ! Actually 4D; first dimension
      ! is path length, same as Eta(:).  Second dimension is Freq X Zeta X Phi
      ! for the state vector.
    logical, intent(inout) :: Do_Calc(:,:)      ! Same shape as Eta_Array
    integer, intent(inout), optional :: NZ(:,:) ! Row subscripts of nonzeros
    integer, intent(inout), optional :: NNZ(:)  ! How many nonzeros in a column
    logical, intent(in), optional :: DerivFlags(:) ! Compute derivatives w.r.t.
                                                ! these state-vector elements.

    integer :: I, J, K
    integer :: N_Path

    n_path = size(eta)

    ! Clean out old nonzeros in Eta_Array and old true values in Do_Calc
    if ( present(nz) ) then
      call clean_out_nonzeros ( eta_array, do_calc, nz, nnz )
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%nf,eta(i)%v(j)%n,eta(i)%v(j)%np], &
                             & three_d_bounds)
          eta_array(i,k) = eta(i)%v(j)%v
          do_calc(i,k) = .true.
          nnz(k) = nnz(k) + 1
          nz(nnz(k),k) = i
        end do
      end do
    else
      eta_array = 0
      do_calc = .false.
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%nf,eta(i)%v(j)%n,eta(i)%v(j)%np], &
                             & three_d_bounds)
          eta_array(i,k) = eta(i)%v(j)%v
          do_calc(i,k) = .true.
        end do
      end do
    end if

    if ( present(derivFlags) ) then
      do i = 1, n_path
        do j = 1, eta(i)%n
          k = element_position([eta(i)%v(j)%n,eta(i)%v(j)%np], three_d_bounds)
          do_calc(i,k) = do_calc(i,k) .and. derivFlags(k)
        end do
      end do
    end if

  end subroutine Get_Eta_Do_Calc_3

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Do_Calc_m

! $Log$
! Revision 2.2  2017/03/02 00:34:30  vsnyder
! Add Clean_Out_NonZeros and Get_Eta_Do_Calc_3
!
! Revision 2.1  2017/02/04 02:19:41  vsnyder
! Initial commit
!

