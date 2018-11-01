! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Tau_M

  use MLSCommon, only: RP

  implicit NONE
  private

  public :: Destroy_Tau, Dump, Dump_Tau, Get_Tau, Black_Out

  type, public :: Tau_T
    real(rp), pointer :: Tau(:,:) => NULL() ! Path X Frequencies
    integer, pointer :: I_Stop(:) => NULL() ! Amount of path to use, size = #Frequencies
  contains
    procedure :: Destroy_Tau
    generic :: Destroy => Destroy_Tau
    procedure :: Dump_Tau
    generic :: Dump => Dump_Tau
  end type Tau_T

  interface Dump; module procedure Dump_Tau; end interface Dump

  ! Use log(underflow) for real(rp), else use log(default real round-off),
  ! to compute black-out point
  logical, private, parameter :: Underflow_Black_Out = .false.

  real(rp), parameter :: black_out = & ! round-off or underflow for exp
!     merge ( -log(huge(1.0_rp)), log(epsilon(1.0)*1.0_rp), Underflow_Black_Out )
    merge ( -log(huge(1.0_rp)), -15.0_rp, Underflow_Black_Out )

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ----------------------------------------------------  Get_Tau  -----
  subroutine Get_Tau ( Frq_i, gl_inds, more_inds, i_start, e_rflty,  &
                     & del_zeta, alpha_path, ref_cor, incoptdepth, &
                     & tan_pt, ds_dz_gw, tau )

  !{ Evaluate {\tt Tau(j)} = $\tau(s_j,s_m) = \int_{s_j}^{s_m}
  !  \alpha(s) \, \text{d} s$, where $s_m$ is the location of the instrument.
  !  Here, {\tt j = 1} is the instrument location.

  ! Update IncOptDepth with its GL corrections.  Multiply it by the
  ! refractive correction Ref_Cor.  Compute Tau = exp(indefinite sum of
  ! IncOptDepth).  The tangent point is a zero-thickness layer that
  ! doesn't have any IncOptDepth.  Instead, multiply tau(tan_pt) by
  ! e_rflty to get tau(tan_pt+1).

  ! This is just Rad_Tran without Inc_Rad_Path and Radiance calculations

    use GL_Update_Incoptdepth_m
    use MLSCommon, only: RP, IP

  ! inputs

    integer(ip), intent(in) :: Frq_I         ! Which frequency slice in Tau_T?
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    integer(ip), intent(in) :: more_inds(:)  ! Places in the coarse path
                                             ! where GL is needed
    integer, intent(in) :: I_Start           ! Where in path to start integrating
    real(rp), intent(in) :: e_rflty          ! Earth reflectivity value (0--1).
    real(rp), intent(in) :: del_zeta(:)      ! Path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just the part up
                                             ! to the black-out
    real(rp), intent(in) :: alpha_path(:)    ! Absorption coefficient on
                                             ! composite coarse & fine grid.
    real(rp), intent(in) :: ref_cor(:)       ! Refracted to unrefracted path
                                             ! length ratios.
    real(rp), intent(inout) :: incoptdepth(:) ! Incremental path opacities
                                             ! from one-sided layer calculation
                                             ! on output. It is the full
                                             ! integrated layer opacity.
    integer, intent(in) :: tan_pt            ! Tangent point index in IncOptDepth
    real(rp), intent(in) :: ds_dz_gw(:)      ! Path length wrt zeta derivative * gw.

  ! outputs

    type(tau_t), intent(inout) :: tau        ! Transmission function.  inout so
                                             ! as not to clobber the association
                                             ! status of its components. 
                                             ! Initial values not used.

  ! Internals

    integer :: I_Stop, N_Path
    real(rp) :: Total_Opacity

  ! Begin code

  ! see if anything needs to be gl-d

    call GL_update_incoptdepth ( gl_inds, more_inds, i_start, del_zeta, &
                               & alpha_path, ds_dz_gw, ref_cor, &
                               & incoptdepth )

  ! Compute Tau = exp(-indefinite sum of IncOptDepth)

    n_path = size(incoptdepth)
    total_opacity = 0.0_rp

    tau%tau(1:i_start,frq_i) = 1.0_rp

    !{ Compute $\tau_i$ for {\tt i_start} $ < i \leq t$, where $t$ is given by
    !  tan\_pt.
    !  $\tau_i = \exp \left \{ - \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1}
    !  \right \}$.
    !  $\Delta \delta_{j \rightarrow j-1}$ is given by incoptdepth and
    !  $- \sum_{j=2}^i \Delta \delta_{j \rightarrow j-1}$ is given by
    !  total\_opacity.

    do i_stop = i_start+1, min(tan_pt,n_path)
      total_opacity = total_opacity - incoptdepth(i_stop)
      tau%tau(i_stop,frq_i) = exp(total_opacity)
      if ( total_opacity < black_out ) go to 19 ! Won't add anything
    end do
    if ( tan_pt >= n_path ) then
      i_stop = i_stop - 1
      go to 19
    end if

    !{ Account for earth reflectivity at the tangent to the surface:
    !  $\tau_{t + 1} = \Upsilon \tau_t$.

    if ( i_start < tan_pt ) then
      tau%tau(tan_pt+1,frq_i) = e_rflty * tau%tau(tan_pt,frq_i)
      total_opacity = total_opacity + log(e_rflty)
    else
      i_stop = i_start
    end if

    !{ Compute $\tau_i$ for $i > 2 N - t + 1$, where $t$ is given by tan\_pt.\\
    !  $\tau_i = \tau_{2N - t + 1} \exp \left \{ - \sum_{j=2N - t + 1}^i
    !    \Delta \delta_{j-1 \rightarrow j} \right \}$.

    ! i_stop is tan_pt + 1 here.

    if ( total_opacity >= black_out ) then
      do while ( i_stop < n_path )
        total_opacity = total_opacity - incoptdepth(i_stop)
        if ( total_opacity < black_out ) exit ! Won't add anything
        i_stop = i_stop + 1
        tau%tau(i_stop,frq_i) = exp(total_opacity)
      end do
    end if

19  continue
    tau%tau(i_stop+1:n_path,frq_i) = 0.0_rp

    tau%i_stop(frq_i) = i_stop

  end subroutine Get_Tau

! --------------------------------------------------  Destroy_Tau  -----
  subroutine Destroy_Tau ( Tau, What, Where )
    use Allocate_Deallocate, only: Deallocate_test
    class(tau_t), intent(inout) :: Tau
    character(len=*), intent(in) :: What, Where
    call deallocate_test ( tau%tau, what//"%Tau", where )
    call deallocate_test ( tau%i_stop, what//"%I_Stop", where )
  end subroutine Destroy_Tau

! -----------------------------------------------------  Dump_Tau  -----
  subroutine Dump_Tau ( Tau, NoFreqs, What, I_Start )
    use Dump_0, only: Dump
    use Output_m, only: Output

    class(tau_t), intent(in) :: Tau
    integer, intent(in) :: NoFreqs ! Number of frequences
    character(len=*), intent(in), optional :: What
    integer, intent(in), optional :: I_Start

    integer :: I, My_Start

    if ( present(what) ) call output ( what, advance='yes' )
    my_Start = 1
    if ( present(i_start) ) my_Start = i_start

    do i = 1, noFreqs
      call dump ( tau%tau(my_start:tau%i_stop(i),i) )
    end do

  end subroutine Dump_Tau

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Tau_M

! $Log$
! Revision 2.17  2018/11/01 00:47:59  vsnyder
! Add Underflow_Black_Out parameter.  Define Black_Out parameter using
! Merge to select log(Underflow) or -15 as black-out value.
!
! Revision 2.16  2018/05/14 23:29:41  vsnyder
! Make some procedures type bound, use GL_Update_Incoptdepth_m
!
! Revision 2.15  2014/07/18 23:16:08  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.14  2011/09/30 01:53:43  vsnyder
! Correct a comment
!
! Revision 2.13  2010/12/22 21:46:25  vsnyder
! TeXnicalities
!
! Revision 2.12  2010/02/02 01:36:03  vsnyder
! Don't reference incoptdepth and ref_cor before i_start
!
! Revision 2.11  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.10  2009/06/13 01:13:02  vsnyder
! Specify start and end of path
!
! Revision 2.9  2007/12/04 23:39:19  vsnyder
! Make black_out public
!
! Revision 2.8  2007/11/07 02:34:32  vsnyder
! Correct a comment
!
! Revision 2.7  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.6  2005/11/01 23:02:21  vsnyder
! PFA Derivatives
!
! Revision 2.5  2005/10/24 20:14:41  vsnyder
! Tau after black out is zero, not one
!
! Revision 2.4  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2005/03/28 20:23:41  vsnyder
! Add nFreqs argument to dump routine
!
! Revision 2.2  2005/03/03 02:08:32  vsnyder
! Add Destroy, Dump routines
!
! Revision 2.1  2004/10/07 17:34:15  vsnyder
! Initial commit
!
