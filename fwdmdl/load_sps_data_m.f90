! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LOAD_SPS_DATA_M
  use MLSCommon, only: RP, IP
  use Units, only: Deg2Rad
  use Intrinsic, only: L_VMR
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use Molecules, only: spec_tags
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  implicit none

  Private
  Public :: load_sps_data
!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
 "$Id$"
  character (LEN=*), parameter :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-------------------------------------------------------------------

  subroutine load_sps_data(fwdModelIn,fwdModelExtra,molecules,        &
    &   f_len,n_f_phi,n_f_zeta,h2o_ind,       &
    &   no_z,no_phi,lin_log,z_basis,phi_basis,sps_values,             &
    &   n_dw_phi,n_dw_zeta,n_dn_phi,n_dn_zeta,        &
    &   n_dv_phi,n_dv_zeta,no_z_dw,no_phi_dw,no_z_dn,no_phi_dn, &
    &   no_z_dv,no_phi_dv,z_basis_dw,phi_basis_dw,z_basis_dn,         &
    &   phi_basis_dn,z_basis_dv,phi_basis_dv)

    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    integer, intent(in)    :: molecules(:)

    integer, intent(out) :: f_len
    integer, intent(out) :: h2o_ind
    integer, intent(out) :: n_f_phi
    integer, intent(out) :: n_f_zeta

    integer, pointer :: no_z(:), no_phi(:)

    logical, pointer :: lin_log(:)

    real(rp), pointer :: z_basis(:), phi_basis(:), sps_values(:)

    character(LEN=3), parameter :: WNV='+++'

!   character(LEN=*), optional, intent(in) :: WNV
!   type(Spect_der_T), optional, intent(in) :: Spect_Der(:)

    integer, intent(out) :: n_dw_phi, n_dw_zeta
    integer, intent(out) :: n_dn_phi, n_dn_zeta
    integer, intent(out) :: n_dv_phi, n_dv_zeta

    integer, pointer :: no_z_dw(:), no_phi_dw(:)
    integer, pointer :: no_z_dn(:), no_phi_dn(:)
    integer, pointer :: no_z_dv(:), no_phi_dv(:)

    real(rp), pointer :: z_basis_dw(:), phi_basis_dw(:)
    real(rp), pointer :: z_basis_dn(:), phi_basis_dn(:)
    real(rp), pointer :: z_basis_dv(:), phi_basis_dv(:)

    ! Local variables:

    Character(len=1) :: CA
    integer :: i,j,k,l,m,n,kz,kp,Spectag,n_sps,accum_z_dw,accum_p_dw, &
           &   accum_z_dn,accum_p_dn,accum_z_dv,accum_p_dv,j_dw,l_dw,&
           &   j_dn,l_dn,j_dv,l_dv

    type (VectorValue_T), pointer :: F             ! An arbitrary species
!
    !******************* LOAD SPECIES DATA ************

    n_sps = size ( molecules )
    call Allocate_test ( no_z, n_sps, 'no_z', ModuleName )
    call Allocate_test ( no_phi, n_sps, 'no_phi', ModuleName )
    call Allocate_test ( lin_log, n_sps, 'lin_log', ModuleName )

    h2o_ind = 0
    do i = 1 , n_sps
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=molecules(i) )
      no_phi(i) = f%template%noInstances
      no_z(i) = f%template%noSurfs
      if(spec_tags(molecules(i)) == 18003) h2o_ind = i
    end do
!
    f_len = sum(no_z * no_phi)
    n_f_phi = sum(no_phi)
    n_f_zeta = sum(no_z)
    call Allocate_test ( sps_values, f_len, 'sps_values', ModuleName )
    call Allocate_test ( z_basis, n_f_zeta, 'z_basis', ModuleName )
    call Allocate_test ( phi_basis, n_f_phi, 'phi_basis', ModuleName )
!
    j = 1
    l = 1
    f_len = 1
    do i = 1 , n_sps
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=molecules(i) )
      k = j + no_phi(i)
      n = l + no_z(i)
      m = f_len + no_z(i) * no_phi(i)
      z_basis(l:n-1) = f%template%surfs(:,1)
      phi_basis(j:k-1) = f%template%phi(1,:) * Deg2Rad
      sps_values(f_len:m-1) = reshape(f%values(1:no_z(i),1:no_phi(i)), &
        &  (/no_z(i)*no_phi(i)/))
      if (f%template%logBasis) then
        lin_log(i) = .true.
        sps_values(f_len:m-1) = log(sps_values(f_len:m-1))
      else
        lin_log(i) = .false.
      endif
      j = k
      l = n
      f_len = m
    enddo
    f_len = f_len - 1
!
    !******************* LOAD SPECTRAL SPECIES DATA ****************
!
    !*** if (.not. associated(Spect_der) ) return
    return
!
    if(index(WNV,'W') > 0) then
      call Allocate_test ( no_z_dw, n_sps, 'no_z_dw', ModuleName )
      call Allocate_test ( no_phi_dw, n_sps, 'no_phi_dw', ModuleName )
      no_z_dw = 0
      no_phi_dw = 0
    endif
!
    if(index(WNV,'N') > 0) then
      call Allocate_test ( no_z_dn, n_sps, 'no_z_dn', ModuleName )
      call Allocate_test ( no_phi_dn, n_sps, 'no_phi_dn', ModuleName )
      no_z_dn = 0
      no_phi_dn = 0
    endif
!
    if(index(WNV,'V') > 0) then
      call Allocate_test ( no_z_dv, n_sps, 'no_z_dv', ModuleName )
      call Allocate_test ( no_phi_dv, n_sps, 'no_phi_dv', ModuleName )
      no_z_dv = 0
      no_phi_dv = 0
    endif
!
    accum_z_dw = 0
    accum_p_dw = 0
    accum_z_dn = 0
    accum_p_dn = 0
    accum_z_dv = 0
    accum_p_dv = 0
!
    do i = 1, n_sps
      m = 0
      Spectag = spec_tags(molecules(i))
      do
        m = m + 1
        !  *** if(Spect_Der(m)%Spectag == Spectag) exit
        exit
        if(m == 3*n_sps) exit
      end do
      !  *** if(Spect_Der(m)%Spectag /= Spectag) cycle
      CA = '@' ! *** Spect_Der(m)%type
      kz = 0 ! Spect_Der(m)%no_zeta_values
      kp = 0 ! Spect_Der(m)%no_phi_values
      select case ( CA )
        case ( 'W' )
          no_z_dw(i) = kz
          no_phi_dw(i) = kp
          accum_z_dw = accum_z_dw + kz
          accum_p_dw = accum_p_dw + kp
        case ( 'N' )
          no_z_dn(i) = kz
          no_phi_dn(i) = kp
          accum_z_dn = accum_z_dn + kz
          accum_p_dn = accum_p_dn + kp
        case ( 'V' )
          no_z_dv(i) = kz
          no_phi_dv(i) = kp
          accum_z_dv = accum_z_dv + kz
          accum_p_dv = accum_p_dv + kp
      end select
    end do
!
    if(accum_z_dw * accum_p_dw < 1) then
      call Deallocate_test ( no_z_dw, 'no_z_dw', ModuleName )
      call Deallocate_test ( no_phi_dw, 'no_phi_dw' ,ModuleName )
    else
      n_dw_phi = accum_p_dw
      n_dw_zeta = accum_z_dw
      call Allocate_test(z_basis_dw,n_dw_zeta,'z_basis_dw',ModuleName)
      call Allocate_test(phi_basis_dw,n_dw_phi,'phi_basis_dw',ModuleName)
    endif
!
    if(accum_z_dn * accum_p_dn < 1) then
      call Deallocate_test ( no_z_dn, 'no_z_dn', ModuleName )
      call Deallocate_test ( no_phi_dn, 'no_phi_dn' ,ModuleName )
    else
      n_dn_phi = accum_p_dn
      n_dn_zeta = accum_z_dn
      call Allocate_test(z_basis_dn,n_dn_zeta,'z_basis_dn',ModuleName)
      call Allocate_test(phi_basis_dn,n_dn_phi,'phi_basis_dn',ModuleName)
    endif
!
    if(accum_z_dv * accum_p_dv < 1) then
      call Deallocate_test ( no_z_dv, 'no_z_dv', ModuleName )
      call Deallocate_test ( no_phi_dv, 'no_phi_dv' ,ModuleName )
    else
      n_dv_phi = accum_p_dv
      n_dv_zeta = accum_z_dv
      call Allocate_test(z_basis_dv,n_dv_zeta,'z_basis_dv',ModuleName)
      call Allocate_test(phi_basis_dv,n_dv_phi,'phi_basis_dv',ModuleName)
    endif
!
    j_dw = 1
    l_dw = 1
    j_dn = 1
    l_dn = 1
    j_dv = 1
    l_dv = 1
!
    do i = 1, n_sps
      m = 0
      Spectag = spec_tags(molecules(i))
      do
        m = m + 1
        !  *** if(Spect_Der(m)%Spectag == Spectag) exit
        exit
        if(m == 3*n_sps) exit
      end do
      !  *** if(Spect_Der(m)%Spectag /= Spectag) cycle
      CA = '@' ! *** Spect_Der(m)%type
      select case ( CA )
        case ( 'W' )
          kz = no_z_dw(i)
          kp = no_phi_dw(i)
          if(kp*kz > 0) then
            k = j_dw + kp
            n = l_dw + kz
            !** z_basis_dw(l_dw:n-1)=Spect_Der(m)%zeta_basis(1:kz)
            !** phi_basis_dw(j_dw:k-1)=Spect_Der(m)%phi_basis(1:kp)*Deg2Rad
            j_dw = k
            l_dw = n
          endif
        case ( 'N' )
          kz = no_z_dn(i)
          kp = no_phi_dn(i)
          if(kp*kz > 0) then
            k = j_dn + kp
            n = l_dn + kz
            !** z_basis_dn(l_dn:n-1)=Spect_Der(m)%zeta_basis(1:kz)
            !** phi_basis_dn(j_dn:k-1)=Spect_Der(m)%phi_basis(1:kp)*Deg2Rad
            j_dn = k
            l_dn = n
          endif
        case ( 'V' )
          kz = no_z_dv(i)
          kp = no_phi_dv(i)
          if(kp*kz > 0) then
            k = j_dv + kp
            n = l_dv + kz
            !** z_basis_dv(l_dv:n-1)=Spect_Der(m)%zeta_basis(1:kz)
            !** phi_basis_dv(j_dv:k-1)=Spect_Der(m)%phi_basis(1:kp)*Deg2Rad
            j_dv = k
            l_dv = n
          endif
      end select
    end do
!
  end subroutine load_sps_data

end module LOAD_SPS_DATA_M
! $Log$
! Revision 1.1.2.10  2001/09/12 21:38:51  zvi
! Added CVS stuff
!
! Revision 1.1.2.9  2001/09/08 23:11:41  zvi
! Bug fixed..
!
! Revision 1.1.2.8  2001/09/08 22:43:37  zvi
! Bug fixed..
!
! Revision 1.1.2.7  2001/09/08 22:37:11  zvi
! Eliminatin molecule coeff. stuff
!
! Revision 1.1.2.6  2001/09/08 20:22:00  zvi
! Fixing a bug..
!
! Revision 1.1.2.5  2001/09/08 20:19:48  zvi
! Developing code..
!
! Revision 1.1.2.4  2001/09/08 18:30:49  zvi
! Working on developement version
!
! Revision 1.1.2.3  2001/09/07 20:41:42  livesey
! Messing around.
!
! Revision 1.1.2.2  2001/09/07 20:16:28  livesey
! Changed stuff to lower case
!
! Revision 1.1.2.1  2001/09/07 19:56:41  zvi
! New module..
!
! Revision 1.0  2001/09/06 13:07:09  zvi
