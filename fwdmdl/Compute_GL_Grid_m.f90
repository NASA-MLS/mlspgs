! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Compute_GL_Grid_M

  implicit NONE
  private
  public :: Compute_GL_Grid

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------  Compute_GL_Grid  -----

  subroutine Compute_GL_Grid ( FwdModelConf, Temp, Qtys, &
    &                          Nlvl, MaxVert, P_GLgrid, Z_GLgrid, &
    &                          Tan_Inds, Tan_Press )

  ! Compute the pressure and zeta GL grids.  Compute Tan_Inds and Tan_Press
  ! because they depend on an intermediate result of the GL grid calculation
  ! (rec_tan_inds), that isn't needed anywhere else.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelVectorTools, only: QtyStuff_T
    use GLnp, only: NG, GX
    use Make_Z_Grid_M, only: Make_Z_Grid
    use MLSCommon, only: RP, IP
    use VectorsModule, only: VectorValue_T

  ! Inputs:
    type (forwardModelConfig_T), intent(in) :: fwdModelConf
    type (vectorValue_T), intent(in) :: TEMP      ! Temperature component of state vector
    type (qtyStuff_t), intent(in) :: Qtys(:)      ! Array of pointers to Qty's.

  ! Outputs
    integer, intent(out) :: Nlvl                  ! Levels in coarse grid
    integer, intent(out) :: MaxVert               ! Levels in find grid

  ! Would be intent(out) if they weren't pointers.  First thing here
  ! is to nullify them.
    real(rp), dimension(:), pointer :: P_GLgrid   ! Pressure on glGrid surfs
    real(rp), dimension(:), pointer :: Z_GLgrid   ! Zeta on glGrid surfs
    integer, dimension(:), pointer :: Tan_Inds    ! Index of tangent grid into gl grid
    real(rp), dimension(:), pointer :: Tan_Press

  ! Local variables
    integer :: J
    integer, parameter :: Ngp1 = Ng+1  ! NG + 1
    integer :: NLM1                               ! NLVL - 1
    integer :: No_Tan_Hts
    integer(ip), dimension(:), pointer :: Rec_Tan_Inds ! recommended tangent
      !                                             point indices from make_z_grid
    integer :: SPS_I
    integer :: Z_All_Prev, Z_All_Size
    real(rp), dimension(:), pointer :: Z_all  ! mass storage of representation
      !                                  bases for z_grid determination
    real(rp), dimension(:), pointer :: Z_psig ! recommended PSIG for
      !                                  radiative transfer calculations
      ! THIS VARIABLE REPLACES FwdModelConf%integrationGrid%surfs

    nullify ( p_glgrid, rec_tan_inds, tan_inds, tan_press, z_all, z_glgrid, &
      & z_psig )

! Insert automatic preselected integration gridder here. Need to make a
! large concatenated vector of bases and pointings.

! Calculate size of z_all and allocate it

    z_all_size = temp%template%nosurfs+2 + &
      & Size(FwdModelConf%integrationGrid%surfs)
    if ( associated(FwdModelConf%tangentGrid) ) &
      & z_all_size = z_all_size + FwdModelConf%tangentGrid%nosurfs
    do sps_i = 1 , size(qtys)
      z_all_size = z_all_size + qtys(sps_i)%qty%template%nosurfs
    end do
    call allocate_test ( z_all, z_all_size, 'z_all', moduleName )

! Fill in z_all
! the -3.000 is a designated "surface" value

    z_all_prev = temp%template%nosurfs+2
    z_all(1) = -3.000_rp
    z_all(2:z_all_prev-1) = temp%template%surfs(:,1)
    z_all(z_all_prev) = 4.000_rp

    ! Add the original Integration Grid:
    z_all_size = z_all_prev + Size(FwdModelConf%integrationGrid%surfs)
    z_all(z_all_prev+1:z_all_size) = FwdModelConf%integrationGrid%surfs(:,1)
    z_all_prev = z_all_size

    if ( associated(FwdModelConf%tangentGrid) ) then
      ! if pointing grid is associated concatenate it to the state vector
      z_all_size = z_all_prev + FwdModelConf%tangentGrid%nosurfs
      z_all(z_all_prev+1:z_all_size) = FwdModelConf%tangentGrid%surfs(:,1)
      z_all_prev = z_all_size
    end if

    do sps_i = 1 , size(qtys)
      z_all_size = z_all_size + qtys(sps_i)%qty%template%nosurfs
      z_all(z_all_prev+1:z_all_size) = qtys(sps_i)%qty%template%surfs(:,1)
      z_all_prev = z_all_size
    end do

! Now, create the final grid and discard the temporary array:

    call make_z_grid ( z_all, z_psig, rec_tan_inds )
    call deallocate_test ( z_all, 'z_all', moduleName )

! note that z_psig(1) is the designated surface
    Nlvl = SIZE(z_psig)
    NLm1 = Nlvl - 1
    maxVert = NLm1 * Ngp1 + 1

! Allocate GL grid stuff

    call allocate_test ( z_glGrid,    maxVert, 'z_glGrid', moduleName )
    call allocate_test ( p_glGrid,    maxVert, 'p_glGrid', moduleName )

! From the selected integration grid pressures define the GL pressure grid:

    z_glgrid(1:maxVert-1) = reshape ( &
      ! Midpoint of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) + z_psig(1:Nlm1)),1,Ngp1) + &
      ! Half length of integration grid intervals:
      & spread(0.5_rp * (z_psig(2:Nlvl) - z_psig(1:Nlm1)),1,Ngp1) * &
      ! Gauss points (with -1 at front):
      & spread((/-1.0_rp,Gx(1:Ng)/),2,NLm1), (/maxVert-1/))
    z_glgrid(maxVert) = z_psig(Nlvl)
    p_glgrid = 10.0_rp**(-z_glgrid)

    call deallocate_test ( z_psig, 'z_psig', moduleName )

! Allocate tan_inds and tan_press.

    j = COUNT(fwdModelConf%tangentGrid%surfs < (z_glgrid(1) - 0.0001_rp))
    no_tan_hts = Nlvl + j

    call allocate_test ( tan_inds,  no_tan_hts, 'tan_inds',  moduleName )
    call allocate_test ( tan_press, no_tan_hts, 'tan_press', moduleName )

! Compute tan_inds from rec_tan_inds

    tan_inds(1:j) = 1
    tan_inds(j+1:no_tan_hts) = (rec_tan_inds - 1) * Ngp1 + 1

    call deallocate_test ( rec_tan_inds, 'rec_tan_inds', moduleName )

! Compute tan_press from fwdModelConf%tangentGrid%surfs and z_glgrid

    tan_press(1:j) = fwdModelConf%tangentGrid%surfs(1:j,1)
    tan_press(j+1:no_tan_hts) = z_glgrid(tan_inds(j+1:no_tan_hts))

  end subroutine Compute_GL_Grid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Compute_GL_Grid_M

! $Log$
! Revision 2.3  2003/06/20 19:35:17  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.2  2003/05/20 00:06:23  vsnyder
! Remove stuff not used by FullForwardModel
!
! Revision 2.1  2003/05/15 20:25:23  vsnyder
! Initial commit
!
