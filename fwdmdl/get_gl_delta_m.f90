! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_GL_DELTA_M

  implicit NONE
  private
  public :: GET_GL_DELTA
  
  interface Get_gl_delta
    module procedure Get_gl_delta_scalar, Get_gl_delta_polarized
  end interface  

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

! ------------------------------------------  Get_gl_delta_scalar  -----

  subroutine Get_gl_delta_scalar ( gl_inds, do_gl, z_path_c, &
                        & alpha_path_c, alpha_path_gl, ds_dh_gl, dh_dz_gl,&
                        & gl_delta )

  ! Compute the correction to the rectangle rule resulting from Gauss-
  ! Legendre quadrature on the fine grid.

    use DO_DELTA_M, ONLY: PATH_OPACITY
    use GLNP, ONLY: Ng
    use MLSCommon, only: RP, IP

    integer(ip), intent(in) :: Gl_inds(:)   ! Gauss-Legendre grid indices
    logical, intent(in) :: Do_gl(:)         ! path flag indicating where to do
  !                                           GL integrations.
    real(rp), intent(in) :: Z_path_c(:)     ! path -log(P) on coarse grid.
    real(rp), intent(in) :: Alpha_path_c(:) ! absorption coefficient
    real(rp), intent(in) :: Alpha_path_gl(:) ! absorption coefficient on GL
  !                                           grid.
    real(rp), intent(in) :: Ds_dh_gl(:)     ! path length wrt height derivative
  !                                           on GL grid.
    real(rp), intent(in) :: Dh_dz_gl(:)     ! path height wrt zeta derivative
  !                                           on GL grid.
  ! outputs
    real(rp), intent(out) :: Gl_delta(:)    ! gl corrections to selected slabs

  ! local variables
    real(rp) :: del_zeta( size(gl_inds)/NG )
    logical :: Do_It                    ! Saw one do_gl true
    integer(ip) :: more_inds( size(gl_inds)/NG )
    integer(ip) :: i, n_path, p_i

  ! see if anything needs to be gl-d

    do_it = .false.
    n_path = size(z_path_c)
    i = 1
    do p_i = 1, n_path
      if ( do_gl(p_i) ) then
        do_it = .true.
        more_inds(i) = p_i
        if ( p_i > n_path/2 ) then
          del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
        else
          del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
        end if
        i = i + 1
      end if
    end do

    if ( do_it ) &
      & call path_opacity ( del_zeta,  &
               &  alpha_path_c(more_inds), &
               &  alpha_path_gl, ds_dh_gl(gl_inds),  &
               &  dh_dz_gl(gl_inds), gl_delta )

  end subroutine Get_gl_delta_scalar

! ---------------------------------------  Get_gl_delta_polarized  -----

  subroutine Get_gl_delta_polarized ( indices_c, gl_inds, do_gl, z_path, &
                        & alpha_path_c, alpha_path_gl,  ds_dh_gl, dh_dz_gl,&
                        & gl_delta)

    use DO_DELTA_M, ONLY: POLARIZED_PATH_OPACITY
    use GLNP, ONLY: Ng
    use MLSCommon, only: RP, IP

    integer(ip), intent(in) :: Indices_c(:) ! coarse grid indices
    integer(ip), intent(in) :: Gl_inds(:)    ! Gauss-Legendre grid indices
    logical, intent(in) :: Do_gl(:) ! path flag indicating where to do
     !                                gl integrations.
    real(rp), intent(in) :: Z_path(:) ! path -log(P) on input grid.
    complex(rp), intent(in) :: Alpha_path_c(:,:) ! absorption coefficient
    complex(rp), intent(in) :: Alpha_path_gl(:,:) ! absorption coefficient on
     !                                    gl grid.
    real(rp), intent(in) :: Ds_dh_gl(:) ! derivative of path length wrt height
     !                                    on gl grid.
    real(rp), intent(in) :: Dh_dz_gl(:) ! derivative of path height wrt zeta on
     !                                    gl grid.
  ! outputs
    complex(rp), intent(out) :: Gl_delta(:,:) ! gl corrections to selected panels
 
  ! local variables
    real(rp) :: del_zeta( size(alpha_path_c,2) ) 
    integer(ip) :: i, n_path, p_i

  ! see if anything needs to be gl-d

    n_path = size(indices_c)
    i = 0
    do p_i = 1, n_path
      if ( do_gl(p_i) ) then
        i = i + 1
        if ( p_i > n_path/2 ) then
          del_zeta(p_i) = z_path(indices_c(p_i+1)) - z_path(indices_c(p_i))
        else
          del_zeta(p_i) = z_path(indices_c(p_i-1)) - z_path(indices_c(p_i))
        end if
      else
        del_zeta(p_i) = 0.0_rp
        gl_delta(:,p_i) = 0.0_rp
      end if
    end do
    if ( i > 0 ) &
      & call polarized_path_opacity ( del_zeta, do_gl, &
               &  alpha_path_c, &
               &  alpha_path_gl, ds_dh_gl(gl_inds),  &
               &  dh_dz_gl(gl_inds), gl_delta )

  end subroutine Get_gl_delta_polarized

! --------------------------------------------------  not_used_here  -----
  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module GET_GL_DELTA_M

!$Log$
!Revision 2.4  2003/05/20 00:05:13  vsnyder
!Cosmetic changes
!
!Revision 2.3  2003/05/05 23:00:25  livesey
!Merged in feb03 newfwm branch
!
!Revision 2.2.2.2  2003/03/05 03:40:35  vsnyder
!More polarized work
!
!Revision 2.2.2.1  2003/03/01 02:40:21  vsnyder
!Cosmetic changes
!
!Revision 2.2  2003/02/07 00:22:35  michael
!it will compile now
!
!Revision 2.1  2003/02/06 23:24:33  michael
!initial commit
!
