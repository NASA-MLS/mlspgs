! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Do_Calc_Indexed_m

  implicit NONE
  private
  public :: Get_Do_Calc_Indexed
  public :: Get_Do_Calc_Indexed_Coarse, Get_Do_Calc_Indexed_GL

  interface Get_Do_Calc_Indexed
    module procedure Get_Do_Calc_Indexed_Coarse, Get_Do_Calc_Indexed_GL
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ---------------------------------  Get_Do_Calc_Indexed_Coarse  -----
  subroutine Get_Do_Calc_Indexed_Coarse  ( N, Tan_Pt_C, Do_Calc_All, &
    & N_Inds, Inds, debug )

  ! Store L in Inds if Do_Calc_All is set for either boundary of layer L.
  ! Before the tangent, layer L is from boundary L-1 to L.
  ! After the tangent point, layer L is from boundary L to L+1.
  
    use GLNP, ONLY: Ng, NGP1

    integer, intent(in) :: N              ! sizes on coarse grid
    integer, intent(in) :: Tan_Pt_C       ! Index of tangent point in coarse grid

    logical, contiguous, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer,             intent(out) :: N_Inds        ! How many Inds produced
    integer, contiguous, intent(out) :: Inds(:,:)     ! Nx2, Inds(:,1) is the
                                          ! layer index and the index of the
                                          ! boundary nearest the tangent point;
                                          ! Inds(:,2) is the other boundary index.

    logical, optional, intent(in) :: Debug
    integer :: M, P_I, P1, P2
    logical :: DebugHere

    DebugHere = .false.
    if ( present(Debug) ) DebugHere = Debug
    n_inds = 0
    p1 = 2
    p2 = tan_pt_c
    do m = -1, 1, 2
      do p_i = p1, p2
        ! Interpolating coefficient at either coarse boundary is nonzero
        if ( do_calc_all(ngp1*p_i-ng) .or. &
             do_calc_all(ngp1*(p_i+m) - ng) ) then
          n_inds = n_inds + 1
          inds(n_inds,1:2) = [p_i,p_i+m*ngp1]
!           if ( DebugHere ) then
!             print *, 'ngp1*p_i-ng, ngp1*(p_i+m) - ng ', ngp1*p_i-ng, ngp1*(p_i+m) - ng
!             print *, 'do_[1], do_[2] ', &
!               & do_calc_all(ngp1*p_i-ng), &
!               & do_calc_all(ngp1*(p_i+m) - ng)
!             print *, 'm, p_i, p_i+m*ngp1, n_inds ', m, p_i, p_i+m*ngp1, n_inds
!           endif
        end if
      end do ! p_i
      p1 = tan_pt_c + 1
      p2 = n - 1
    end do ! m
    
    
    if ( DebugHere ) then
      p1 = 2
      p2 = tan_pt_c
      do m = -1, 1, 2
      if ( DebugHere ) print *, 'p1, p2 ', p1, p2
        do p_i = p1, p2
          print *, 'ngp1*p_i-ng, ngp1*(p_i+m) - ng ', ngp1*p_i-ng, ngp1*(p_i+m) - ng
          print *, 'do_[1], do_[2] ', &
            & do_calc_all(ngp1*p_i-ng), &
            & do_calc_all(ngp1*(p_i+m) - ng)
          print *, 'm, p_i, p_i+m*ngp1, n_inds ', m, p_i, p_i+m*ngp1, n_inds
        enddo
        p1 = tan_pt_c + 1
        p2 = n - 1
      enddo
    endif

  end subroutine Get_Do_Calc_Indexed_Coarse

  ! -------------------------------------  Get_Do_Calc_Indexed_GL  -----
  subroutine Get_Do_Calc_Indexed_GL ( N, Tan_Pt_C, Do_Calc_All, F_Inds, Do_GL, &
    & Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_All(1::ngp1), or Do_GL and any of the corresponding
  ! Do_Calc_All(f_inds) flags are set.
  ! Get_Do_Calc_Indexed determines that Get_d_Delta_df needs to integrate a
  ! panel if Do_Calc_All(1::ngp1) is true, which means the interpolating
  ! coefficient is nonzero, or if Do_GL is true and Do_Calc_All(f_inds) is
  ! true, which means that GL was needed, even if some interpolating
  ! coefficient is zero.

    use GLNP, ONLY: Ng, NGP1

    integer, intent(in) :: N              ! sizes on coarse grid
    integer, intent(in) :: Tan_Pt_C       ! Index of tangent point in coarse grid

    logical, contiguous, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, contiguous, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
                                                      ! GL points on panels needing GL.
    logical, contiguous, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, contiguous, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer,             intent(out) :: N_Inds        ! count(do_calc)
    integer, contiguous, intent(out) :: Inds(:,:)     ! Nx2, Inds(:,1) is the
                                          ! layer index and the index of the
                                          ! boundary nearest the tangent point;
                                          ! Inds(:,2) is the other boundary index.


    integer :: I, K, M, P_I, P1, P2

    i = 1 - Ng
    n_inds = 0
    p1 = 2
    p2 = tan_pt_c
    do m = -1, 1, 2
      do p_i = p1, p2
        ! Interpolating coefficient at either coarse boundary is nonzero
        do_calc(p_i) = do_calc_all(ngp1*p_i-ng) .or. &
                       do_calc_all(ngp1*(p_i+m) - ng)
        if ( do_gl(p_i) ) then
          i = i + Ng
          k = f_inds(i)
          do_calc(p_i) = do_calc(p_i) .or. &
            & do_calc_all(k) .or. do_calc_all(k+1) .or. do_calc_all(k+2)
        end if
        if ( do_calc(p_i) ) then
          n_inds = n_inds + 1
          inds(n_inds,1:2) = [p_i,p_i+m*ngp1]
        end if
      end do ! p_i
      p1 = tan_pt_c + 1
      p2 = n - 1
    end do ! m

  end subroutine Get_Do_Calc_Indexed_GL

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

end module Get_Do_Calc_Indexed_m

! $Log$
! Revision 2.5  2023/06/23 20:44:03  pwagner
! In middle of debugging pol fwdmdl
!
! Revision 2.4  2018/11/19 21:52:36  vsnyder
! Correct confusion between coarse and fine path indexing
!
! Revision 2.3  2018/04/17 22:10:18  vsnyder
! Remove unused declarations
!
! Revision 2.2  2017/08/25 22:49:25  pwagner
! Removed :: that NAG hated
!
! Revision 2.1  2017/08/09 19:57:55  vsnyder
! Initial commit
!
