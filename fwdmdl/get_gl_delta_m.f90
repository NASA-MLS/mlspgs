! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_GL_DELTA_M

  use MLSCommon, only: RP, IP
  use GLNP, ONLY: Ng
  use DO_DELTA_M, ONLY: PATH_OPACITY, POLARIZED_PATH_OPACITY

  implicit NONE
  private
  public :: GET_GL_DELTA
  
  interface get_gl_delta
    module procedure scalar_get_gl_delta, polarized_get_gl_delta
  END INTERFACE  
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!------------------------------------------------------  get_gl_delta  -----


  SUBROUTINE scalar_get_gl_delta ( indices_c, gl_inds, do_gl, z_path, &
                        & alpha_path_c, alpha_path_gl,  ds_dh_gl, dh_dz_gl,&
                        & gl_delta)

    INTEGER(ip), intent(in) :: indices_c(:) ! coarse grid indicies
    INTEGER(ip), INTENT(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    logical, intent(in) :: do_gl(:) ! path flag indicating where to do
  !                                   gl integrations.
    real(rp), intent(in) :: z_path(:) ! path -log(P) on input grid.
    real(rp), intent(in) :: alpha_path_c(:) ! absorption coefficient
    real(rp), intent(in) :: alpha_path_gl(:) ! absorption coefficient on gl
  !                                        grid.
    real(rp), intent(in) :: ds_dh_gl(:) ! path length wrt height derivative on
  !                                        gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
  !                                       gl grid.
  ! outputs
    real(rp), intent(out) :: gl_delta(:) ! gl corrections to selected slabs

  ! local variables
    real(rp) :: del_zeta( size(gl_inds)/NG) 
    integer(ip) :: more_inds( size(gl_inds)/NG) 
    integer(ip) :: i, j, n_path, p_i


 
    n_path = size(indices_c)
    if ( count(do_gl) > 0 ) then

  ! see if anything needs to be gl-d


      i = 1
      j = 1
      do p_i = 1, n_path
        if ( do_gl(p_i) ) then
          more_inds(i) = p_i
          if ( p_i > n_path/2 ) then
            del_zeta(i) = z_path(indices_c(p_i+1)) - z_path(indices_c(p_i))
          else
            del_zeta(i) = z_path(indices_c(p_i-1)) - z_path(indices_c(p_i))
          end if
          i = i + 1
          j = j + Ng
        end if
      end do
     call path_opacity ( del_zeta,  &
                &  alpha_path_c(more_inds), &
                &  alpha_path_gl, ds_dh_gl(gl_inds),  &
                &  dh_dz_gl(gl_inds), gl_delta )
     end if
  end subroutine scalar_get_gl_delta

!------------------------------------------------------  polarized_get_gl_delta  -----


  SUBROUTINE polarized_get_gl_delta ( indices_c, gl_inds, do_gl, z_path, &
                        & alpha_path_c, alpha_path_gl,  ds_dh_gl, dh_dz_gl,&
                        & gl_delta)

    INTEGER(ip), intent(in) :: indices_c(:) ! coarse grid indices
    INTEGER(ip), INTENT(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    logical, intent(in) :: do_gl(:) ! path flag indicating where to do
  !                                   gl integrations.
    real(rp), intent(in) :: z_path(:) ! path -log(P) on input grid.
    complex(rp), intent(in) :: alpha_path_c(:,:) ! absorption coefficient
    complex(rp), intent(in) :: alpha_path_gl(:,:) ! absorption coefficient on gl
  !                                        grid.
    real(rp), intent(in) :: ds_dh_gl(:) ! path length wrt height derivative on
  !                                        gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
  !                                       gl grid.
  ! outputs
    complex(rp), intent(out) :: gl_delta(:,:) ! gl corrections to selected slabs
 
  ! local variables
    real(rp) :: del_zeta( size(gl_inds)/NG) 
    integer(ip) :: more_inds( size(gl_inds)/NG) 
    integer(ip) :: i, j, n_path, p_i


 
    n_path = size(indices_c)
    if ( count(do_gl) > 0 ) then

  ! see if anything needs to be gl-d


      i = 1
      j = 1
      do p_i = 1, n_path
        if ( do_gl(p_i) ) then
          more_inds(i) = p_i
          if ( p_i > n_path/2 ) then
            del_zeta(i) = z_path(indices_c(p_i+1)) - z_path(indices_c(p_i))
          else
            del_zeta(i) = z_path(indices_c(p_i-1)) - z_path(indices_c(p_i))
          end if
          i = i + 1
          j = j + Ng
        end if
      end do
     call polarized_path_opacity ( del_zeta,  &
                &  alpha_path_c(:,more_inds), &
                &  alpha_path_gl, ds_dh_gl(gl_inds),  &
                &  dh_dz_gl(gl_inds), gl_delta )
     end if
  end subroutine polarized_get_gl_delta

! --------------------------------------------------  not_used_here  -----
  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module GET_GL_DELTA_M

!$Log$
!Revision 2.2  2003/02/07 00:22:35  michael
!it will compile now
!
!Revision 2.1  2003/02/06 23:24:33  michael
!initial commit
!
