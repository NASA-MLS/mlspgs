! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LOAD_SPS_DATA_M
  use MLSCommon, only: RP, IP
  use Units, only: Deg2Rad
  use Intrinsic, only: L_VMR
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use Molecules, only: spec_tags
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
  implicit none

  Private
  Public :: load_sps_data

  type, public :: Grids_T              ! Fit all Gridding categories
    integer :: tot_no_frq             ! Total No. of entries in frequency grid
    integer :: tot_no_zet             ! Total No. of entries in zeta grid
    integer :: tot_no_phi             ! Total No. of entries in phi  grid
    integer,  pointer :: no_f(:)      ! No. of entries in frq. grid per sps
    integer,  pointer :: no_z(:)      ! No. of entries in zeta grid per sps
    integer,  pointer :: no_p(:)      ! No. of entries in phi  grid per sps
    real(rp), pointer :: frq_basis(:) ! frq  grid entries (dim: tot_no_frq)
    real(rp), pointer :: zet_basis(:) ! zeta grid entries (dim: tot_no_zet)
    real(rp), pointer :: phi_basis(:) ! phi  grid entries (dim: tot_no_phi)
  end type Grids_T

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
 "$Id$"
  character (LEN=*), parameter :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
contains
!-------------------------------------------------------------------

  subroutine load_sps_data(fwdModelIn, fwdModelExtra, molecules, p_len, &
        &    f_len, h2o_ind, lin_log, sps_values, Grids_f, Grids_dw,    &
        &    Grids_dn, Grids_dv)

    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    integer, intent(in)    :: molecules(:)

    integer, intent(out) :: f_len
    integer, intent(out) :: p_len
    integer, intent(out) :: h2o_ind

    type (Grids_T) :: Grids_f   ! All the coordinates
    type (Grids_T) :: Grids_dw  ! All the spectroscopy(W) coordinates
    type (Grids_T) :: Grids_dn  ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_dv  ! All the spectroscopy(V) coordinates

    logical, pointer :: lin_log(:)

    real(rp), pointer :: sps_values(:)

    character(LEN=3), parameter :: WNV='+++'

!   character(LEN=*), optional, intent(in) :: WNV
!   type(Spect_der_T), optional, intent(in) :: Spect_Der(:)

    ! Local variables:

    Character(len=1) :: CA
    integer :: i,j,k,l,m,n,r,s,kz,kp,kf,Spectag,accum_z_dw,accum_p_dw,j_dw, &
           &   j_dn,l_dn,j_dv,l_dw,l_dv,n_f_phi,n_f_zet,n_f_frq,n_sps,s_dw, &
           &   s_dn,s_dv,accum_z_dn,accum_p_dn,accum_z_dv,accum_p_dv,&
           &   accum_f_dw,accum_f_dn,accum_f_dv

    type (VectorValue_T), pointer :: F             ! An arbitrary species
!
    !******************* LOAD SPECIES DATA ************

    n_sps = size ( molecules )

    allocate ( Grids_f%no_z(n_sps), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%no_z' )
    allocate ( Grids_f%no_p(n_sps), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%no_p' )
    allocate ( Grids_f%no_f(n_sps), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%no_f' )

    call Allocate_test ( lin_log, n_sps, 'lin_log', ModuleName )

    f_len = 0
    p_len = 0
    h2o_ind = 0
    do i = 1 , n_sps
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=molecules(i) )
!     kf = f%template%noChans        ! ** ZEBUG
      kf = 1                         ! ** ZEBUG
      kz = f%template%noSurfs
      kp = f%template%noInstances
      Grids_f%no_z(i) = kz
      Grids_f%no_p(i) = kp
      Grids_f%no_f(i) = kf
      p_len = p_len + kz * kp
      f_len = f_len + kz * kp * kf
      if(spec_tags(molecules(i)) == 18003) h2o_ind = i
    end do

    call Allocate_test ( sps_values, f_len, 'sps_values', ModuleName )
!
    n_f_zet = SUM(Grids_f%no_z)
    n_f_phi = SUM(Grids_f%no_p)
    n_f_frq = SUM(Grids_f%no_f)

    Grids_f%tot_no_zet = n_f_zet
    Grids_f%tot_no_phi = n_f_phi
    Grids_f%tot_no_frq = n_f_frq

    allocate ( Grids_f%zet_basis(n_f_zet), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%zet_basis' )
    allocate ( Grids_f%phi_basis(n_f_phi), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%phi_basis' )
    allocate ( Grids_f%frq_basis(n_f_frq), stat=j )
    if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'Grids_f%frq_basis' )
!
    j = 1
    l = 1
    s = 1
    f_len = 1
    do i = 1 , n_sps
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=molecules(i) )
!     kf = f%template%noChans        ! ** ZEBUG
      kf = 1                         ! ** ZEBUG
      kz = f%template%noSurfs
      kp = f%template%noInstances
      k = j + kp
      n = l + kz
      m = s + kf
      r = f_len + kz * kp * kf
      Grids_f%zet_basis(l:n-1) = f%template%surfs(:,1)
! *** ZEBUG: f%template%%frequencies(:) does not exist yet ..
!     Grids_f%frq_basis(s:m-1) = f%template%%frequencies(:)
      Grids_f%frq_basis(s:m-1) = 0.0
      Grids_f%phi_basis(j:k-1) = f%template%phi(1,:) * Deg2Rad
      sps_values(f_len:r-1)=RESHAPE(f%values(1:kz*kf,1:kp),(/kz*kf*kp/))
      if (f%template%logBasis) then
        lin_log(i) = .true.
        sps_values(f_len:r-1) = log(sps_values(f_len:r-1))
      else
        lin_log(i) = .false.
      endif
      j = k
      l = n
      s = m
      f_len = r
    enddo
    f_len = f_len - 1
!
    !******************* LOAD SPECTRAL SPECIES DATA ****************
!
    !*** if (.not. associated(Spect_der) ) return
    return
!
    if(index(WNV,'W') > 0) then
      allocate ( Grids_dw%no_z(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%no_z' )
      allocate ( Grids_dw%no_p(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%no_p' )
      allocate ( Grids_dw%no_f(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%no_f' )
    endif
!
    if(index(WNV,'N') > 0) then
      allocate ( Grids_dn%no_z(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%no_z' )
      allocate ( Grids_dn%no_p(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%no_p' )
      allocate ( Grids_dn%no_f(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%no_f' )
    endif
!
    if(index(WNV,'V') > 0) then
      allocate ( Grids_dv%no_z(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%no_z' )
      allocate ( Grids_dv%no_p(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%no_p' )
      allocate ( Grids_dv%no_f(n_sps), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%no_f' )
    endif
!
    accum_z_dw = 0
    accum_p_dw = 0
    accum_f_dw = 0
    accum_z_dn = 0
    accum_p_dn = 0
    accum_f_dn = 0
    accum_z_dv = 0
    accum_p_dv = 0
    accum_f_dv = 0
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
      kz = 0   ! Spect_Der(m)%no_zeta_values
      kp = 0   ! Spect_Der(m)%no_phi_values
      kf = 0   ! Spect_Der(m)%no_frq_values
      select case ( CA )
        case ( 'W' )
          Grids_dw%no_z(i) = kz
          Grids_dw%no_p(i) = kp
          Grids_dw%no_f(i) = kf
          accum_z_dw = accum_z_dw + kz
          accum_p_dw = accum_p_dw + kp
          accum_f_dw = accum_f_dw + kf
        case ( 'N' )
          Grids_dn%no_z(i) = kz
          Grids_dn%no_p(i) = kp
          Grids_dn%no_f(i) = kf
          accum_z_dn = accum_z_dn + kz
          accum_p_dn = accum_p_dn + kp
          accum_f_dn = accum_f_dn + kf
        case ( 'V' )
          Grids_dv%no_z(i) = kz
          Grids_dv%no_p(i) = kp
          Grids_dv%no_f(i) = kf
          accum_z_dv = accum_z_dv + kz
          accum_p_dv = accum_p_dv + kp
          accum_f_dv = accum_f_dv + kf
      end select
    end do
!
    if(index(WNV,'W') > 0) then
      Grids_dw%tot_no_zet = accum_z_dw
      Grids_dw%tot_no_phi = accum_p_dw
      Grids_dw%tot_no_frq = accum_f_dw
      allocate ( Grids_dw%zet_basis(accum_z_dw), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%zet_basis' )
      allocate ( Grids_dw%phi_basis(accum_p_dw), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%phi_basis' )
      allocate ( Grids_dw%frq_basis(accum_f_dw), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dw%frq_basis' )
    endif
!
    if(index(WNV,'N') > 0) then
      Grids_dn%tot_no_zet = accum_z_dn
      Grids_dn%tot_no_phi = accum_p_dn
      Grids_dn%tot_no_frq = accum_f_dn
      allocate ( Grids_dn%zet_basis(accum_z_dn), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%zet_basis' )
      allocate ( Grids_dn%phi_basis(accum_p_dn), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%phi_basis' )
      allocate ( Grids_dn%frq_basis(accum_f_dn), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dn%frq_basis' )
    endif
!
    if(index(WNV,'V') > 0) then
      Grids_dv%tot_no_zet = accum_z_dv
      Grids_dv%tot_no_phi = accum_p_dv
      Grids_dv%tot_no_frq = accum_f_dv
      allocate ( Grids_dv%zet_basis(accum_z_dv), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%zet_basis' )
      allocate ( Grids_dv%phi_basis(accum_p_dv), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%phi_basis' )
      allocate ( Grids_dv%frq_basis(accum_f_dv), stat=j )
      if ( j /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//'Grids_dv%frq_basis' )
    endif
!
    j_dw = 1
    l_dw = 1
    s_dw = 1
    j_dn = 1
    l_dn = 1
    s_dn = 1
    j_dv = 1
    l_dv = 1
    s_dv = 1
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
          kz = Grids_dw%no_z(i)
          kp = Grids_dw%no_p(i)
          kf = Grids_dw%no_f(i)
          if(kp*kz*kf > 0) then
            k = j_dw + kp
            n = l_dw + kz
            j = s_dw + kf
      !     Grids_dw%zet_basis(l_dw:n-1) = &
      !                         &  Spect_Der(m)%zeta_basis(1:kz)
      !     Grids_dw%phi_basis(j_dw:k-1) = &
      !                         &  Spect_Der(m)%phi_basis(1:kp)
      !     Grids_dw%frq_basis(s_dw:j-1) = &
      !                         &  Spect_Der(m)%frq_basis(1:kf)
            j_dw = k
            l_dw = n
            s_dw = j
          endif
        case ( 'N' )
          kz = Grids_dn%no_z(i)
          kp = Grids_dn%no_p(i)
          kf = Grids_dn%no_f(i)
          if(kp*kz*kf > 0) then
            k = j_dn + kp
            n = l_dn + kz
            j = s_dn + kf
      !     Grids_dn%zet_basis(l_dn:n-1) = &
      !                         &  Spect_Der(m)%zeta_basis(1:kz)
      !     Grids_dn%phi_basis(j_dn:k-1) = &
      !                         &  Spect_Der(m)%phi_basis(1:kp)
      !     Grids_dn%frq_basis(s_dn:j-1) = &
      !                         &  Spect_Der(m)%frq_basis(1:kf)
            j_dn = k
            l_dn = n
            s_dn = j
          endif
        case ( 'V' )
          kz = Grids_dv%no_z(i)
          kp = Grids_dv%no_p(i)
          kf = Grids_dv%no_f(i)
          if(kp*kz*kf > 0) then
            k = j_dv + kp
            n = l_dv + kz
            j = s_dv + kf
      !     Grids_dv%zet_basis(l_dv:n-1) = &
      !                         &  Spect_Der(m)%zeta_basis(1:kz)
      !     Grids_dv%phi_basis(j_dv:k-1) = &
      !                         &  Spect_Der(m)%phi_basis(1:kp)
      !     Grids_dv%frq_basis(s_dv:j-1) = &
      !                         &  Spect_Der(m)%frq_basis(1:kf)
            j_dv = k
            l_dv = n
            s_dn = j
          endif
      end select
    end do
!
  end subroutine load_sps_data

end module LOAD_SPS_DATA_M
! $Log$
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
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
