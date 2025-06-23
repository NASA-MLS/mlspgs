! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Compute_Z_PSIG_m

  implicit NONE
  private
  public :: Compute_Z_PSIG

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!------------------------------------------------  Compute_Z_PSIG  -----

  subroutine Compute_Z_PSIG ( FwdModelConf, Z_PSIG, Observer )

  ! Compute the preselected integration zeta grid.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use ForwardModelConfig, only: ForwardModelConfig_t, QtyStuff_T
    use Make_Z_Grid_M, only: Make_Z_Grid
    use MLSCommon, only: RP
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_T

  ! Inputs:
    type (forwardModelConfig_T), intent(in) :: fwdModelConf

  ! Outputs:
    real(rp), allocatable, intent(out) :: Z_psig(:) ! recommended PSIG for
      !                                  radiative transfer calculations

  ! Optional inputs:
    type (vectorValue_T), optional, intent(in) :: Observer ! Zetas for observers-in-atmosphere

  ! Local variables:
    integer :: Me = -1             ! String index cache for tracing
    type (qtyStuff_t), pointer :: Qtys(:)         ! Array of pointers to Qty's.
    integer :: SPS_I
    integer :: Z_All_Prev, Z_All_Size
    real(rp), allocatable :: Z_all(:)  ! consolidated storage of representation
      !                              bases for z_grid determination

    call trace_begin ( me, 'Compute_Z_PSIG', &
      & cond=toggle(emit) .and. levels(emit) > 2 ) ! set by -f command-line switch

    qtys => fwdModelConf%beta_group%qty

! Insert automatic preselected integration gridder here. Need to make a
! large concatenated vector of bases and pointings.

! Calculate size of z_all and allocate it

    z_all_size = fwdModelConf%temp%qty%template%nosurfs + 2
    if ( associated(FwdModelConf%integrationGrid) ) &
      & z_all_size = z_all_size + FwdModelConf%integrationGrid%nosurfs
    if ( associated(FwdModelConf%tangentGrid) .and. .not. &
      & associated(FwdModelConf%integrationGrid,FwdModelConf%tangentGrid) ) &
      & z_all_size = z_all_size + FwdModelConf%tangentGrid%nosurfs
    do sps_i = 1 , size(qtys)
      z_all_size = z_all_size + qtys(sps_i)%qty%template%nosurfs
    end do
    if ( present(observer) ) z_all_size = z_all_size + observer%template%nosurfs
    call allocate_test ( z_all, z_all_size, 'z_all', moduleName )

! Fill in z_all
! the -3.000 is a designated "surface" value
! the  4.000 is a designated "top of the atmosphere" value

    z_all_prev = fwdModelConf%temp%qty%template%nosurfs + 2
    z_all(1) = -3.000_rp
    z_all(2:z_all_prev-1) = fwdModelConf%temp%qty%template%surfs(:,1)
    z_all(z_all_prev) = 4.000_rp

    if ( associated(FwdModelConf%integrationGrid) ) then
      ! Add the original Integration Grid
      z_all_size = z_all_prev + FwdModelConf%integrationGrid%nosurfs
      z_all(z_all_prev+1:z_all_size) = FwdModelConf%integrationGrid%surfs(:,1)
      z_all_prev = z_all_size
    end if

    if ( associated(FwdModelConf%tangentGrid) .and. .not. &
      & associated(FwdModelConf%integrationGrid,FwdModelConf%tangentGrid) ) then
      ! if pointing grid is associated and not the same as the integration
      ! grid concatenate it to the state vector
      z_all_size = z_all_prev + FwdModelConf%tangentGrid%nosurfs
      z_all(z_all_prev+1:z_all_size) = FwdModelConf%tangentGrid%surfs(:,1)
      z_all_prev = z_all_size
    end if

    do sps_i = 1, size(qtys)
      z_all_size = z_all_size + qtys(sps_i)%qty%template%nosurfs
      z_all(z_all_prev+1:z_all_size) = qtys(sps_i)%qty%template%surfs(:,1)
      z_all_prev = z_all_size
    end do

    if ( present(observer) ) then
      z_all_size = z_all_prev + observer%template%nosurfs
      z_all(z_all_prev+1:z_all_size) = observer%template%surfs(:,1)
      z_all_prev = z_all_size
    end if

! Now, create the final grid and discard the temporary array:

    call make_z_grid ( z_all, z_psig )
    call deallocate_test ( z_all, 'z_all', moduleName )

    call trace_end ( 'Compute_Z_PSIG', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

  end subroutine Compute_Z_PSIG

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Compute_Z_PSIG_m

! $Log$
! Revision 2.14  2016/09/03 00:30:09  vsnyder
! Make Z_All allocatable instead of a pointer
!
! Revision 2.13  2016/05/02 23:32:31  vsnyder
! Get temperature quantity from FwdModelConf
!
! Revision 2.12  2014/08/06 23:24:51  vsnyder
! Remove USE for Switches, which is not referenced
!
! Revision 2.11  2014/04/22 00:36:54  vsnyder
! Add tracing
!
! Revision 2.10  2013/07/13 00:04:58  vsnyder
! Move computation of tangent pressures to Tangent_Pressures
!
! Revision 2.9  2013/06/12 02:18:56  vsnyder
! Make Z_psig allocatable instead of pointer
!
! Revision 2.8  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2008/08/27 19:56:51  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.6  2008/06/26 00:26:30  vsnyder
! Delete QTYS argument, since the info is in the config argument
!
! Revision 2.5  2008/05/20 00:18:07  vsnyder
! Add Observers dummy argument
!
! Revision 2.4  2008/05/01 01:57:51  vsnyder
! Make sure integrationGrid is associated before asking its size
! Don't add both integrationGrid and tangentGrid if they're the same grid
!
! Revision 2.3  2007/01/17 23:48:43  vsnyder
! Don't look at FwdModelConf%tangentGrid if it's not associated
!
! Revision 2.2  2006/07/12 20:52:57  vsnyder
! Remove declarations for unused variables
!
! Revision 2.1  2006/07/07 17:54:56  vsnyder
! Initial commit
!
