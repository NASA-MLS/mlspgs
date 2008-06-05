! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Read_Mie_m

  ! Read the Mie tables produced by the Mie_Tables program.

  use MLSKinds, only: R8

  implicit NONE
  private
  public :: Destroy_Mie, Dump_Mie, Read_Mie

  ! Coordinates
  real(r8), public, pointer :: F_s(:) => NULL(), IWC_s(:) => NULL()
  real(r8), public, pointer :: T_s(:) => NULL(), Theta_s(:) => NULL()
  ! Betas
  real(r8), public, pointer :: Beta_c_a(:,:,:) => NULL()  ! T X IWC X F
  real(r8), public, pointer :: Beta_c_e(:,:,:) => NULL()  ! T X IWC X F
  real(r8), public, pointer :: Beta_c_s(:,:,:) => NULL()  ! T X IWC X F
  ! Beta derivatives
  real(r8), public, pointer :: dBeta_dIWC_c_a(:,:,:) => NULL() ! T X IWC X F
  real(r8), public, pointer :: dBeta_dT_c_a(:,:,:) => NULL()   ! T X IWC X F
  real(r8), public, pointer :: dBeta_dIWC_c_e(:,:,:) => NULL() ! T X IWC X F
  real(r8), public, pointer :: dBeta_dT_c_e(:,:,:) => NULL()   ! T X IWC X F
  real(r8), public, pointer :: dBeta_dIWC_c_s(:,:,:) => NULL() ! T X IWC X F
  real(r8), public, pointer :: dBeta_dT_c_s(:,:,:) => NULL()   ! T X IWC X F
  ! Phase
  real(r8), public, pointer :: P(:,:,:,:) => NULL()       ! T X IWC X Theta X F
  ! Phase Derivatives
  real(r8), public, pointer :: dP_dIWC(:,:,:,:) => NULL() ! T X IWC X Theta X F
  real(r8), public, pointer :: dP_dT(:,:,:,:) => NULL()   ! T X IWC X Theta X F

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------------  Destroy_Mie  -----
  subroutine Destroy_Mie
  ! Deallocate the Mie tables

    use Allocate_Deallocate, only: E_dp, MEMORY_UNITS, Test_DeAllocate

    integer :: Status
    real :: S

    if ( associated(f_s) ) then
      s = (e_dp * size(f_s) ) / MEMORY_UNITS
      deallocate ( f_s, stat=status )
      call test_deallocate ( status, moduleName, 'F_s', s )
      s = (e_dp * size(iwc_s) ) / MEMORY_UNITS
      deallocate ( iwc_s, stat=status )
      call test_deallocate ( status, moduleName, 'iwc_s', s )
      s = (e_dp * size(t_s) ) / MEMORY_UNITS
      deallocate ( t_s, stat=status )
      call test_deallocate ( status, moduleName, 't_s', s )
      s = (e_dp * size(theta_s) ) / MEMORY_UNITS
      deallocate ( theta_s, stat=status )
      call test_deallocate ( status, moduleName, 'theta_s', s )
      s = (e_dp * size(beta_c_a) ) / MEMORY_UNITS
      deallocate ( beta_c_a, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_a', s )
      s = (e_dp * size(beta_c_e) ) / MEMORY_UNITS
      deallocate ( beta_c_e, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_e', s )
      s = (e_dp * size(beta_c_s) ) / MEMORY_UNITS
      deallocate ( beta_c_s, stat=status )
      call test_deallocate ( status, moduleName, 'beta_c_s', s )
      s = (e_dp * size(P) ) / MEMORY_UNITS
      deallocate ( P, stat=status )
      call test_deallocate ( status, moduleName, 'P', s )
      if ( associated(dP_dIWC) ) then
        s = (e_dp * size(dBeta_dIWC_c_a) ) / MEMORY_UNITS
        deallocate ( dBeta_dIWC_c_a, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_a', s )
        s = (e_dp * size(dBeta_dT_c_a) ) / MEMORY_UNITS
        deallocate ( dBeta_dT_c_a, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_a', s )
        s = (e_dp * size(dBeta_dIWC_c_e) ) / MEMORY_UNITS
        deallocate ( dBeta_dIWC_c_e, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_e', s )
        s = (e_dp * size(dBeta_dT_c_e) ) / MEMORY_UNITS
        deallocate ( dBeta_dT_c_e, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_e', s )
        s = (e_dp * size(dBeta_dIWC_c_s) ) / MEMORY_UNITS
        deallocate ( dBeta_dIWC_c_s, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dIWC_c_s', s )
        s = (e_dp * size(dBeta_dT_c_s) ) / MEMORY_UNITS
        deallocate ( dBeta_dT_c_s, stat=status )
        call test_deallocate ( status, moduleName, 'dBeta_dT_c_s', s )
        s = (e_dp * size(dP_dIWC) ) / MEMORY_UNITS
        deallocate ( dP_dIWC, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dIWC', s )
        s = (e_dp * size(dP_dT) ) / MEMORY_UNITS
        deallocate ( dP_dT, stat=status )
        call test_deallocate ( status, moduleName, 'dP_dT', s )
      end if
    end if

  end subroutine Destroy_Mie

  ! ---------------------------------------------------  Dump_Mie  -----
  subroutine Dump_Mie ( Details )
  ! Dump the Mie tables
  ! Details <= 0 => F_s, T_s, IWC_s, Theta_s
  !          = 1 => 0 + P
  !         >= 2 => 1 + dP_dIWC + dP_dT
    use Dump_0, only: Dump
    use Output_m, only: Output
    use Units, only: Rad2Deg
    integer, intent(in) :: Details
    integer :: I_f, I_theta
    if ( .not. associated(f_s) ) return
    call dump ( f_s, name='F_s (GHz)' )
    call dump ( T_s, name='T_s (K)' )
    call dump ( IWC_s, name='IWC_s (actually log10 IWC)' )
    call dump ( Rad2Deg*theta_s, name='Theta_s (Degrees)' )
    if ( details <= 0 ) return
    do i_f = 1, size(f_s)
      call output ( f_s(i_f), before='Beta_c_a(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_a(:,:,i_f) )
      call output ( f_s(i_f), before='Beta_c_e(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_e(:,:,i_f) )
      call output ( f_s(i_f), before='Beta_c_s(T X IWC), F = ', advance='yes' )
      call dump ( beta_c_s(:,:,i_f) )
      if ( associated(dBeta_dIWC_c_a) .and. details >= 2 ) then
        call output ( f_s(i_f), before='dBeta_dIWC_c_a(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_a(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_a(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_a(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dIWC_c_e(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_e(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_e(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_e(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dIWC_c_s(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dIWC_c_s(:,:,i_f) )
        call output ( f_s(i_f), before='dBeta_dT_c_s(T X IWC), F = ', advance='yes' )
        call dump ( dBeta_dT_c_s(:,:,i_f) )
      end if
      do i_theta = 1, size(theta_s)
        call output ( Rad2Deg*theta_s(i_theta), before='P(T X IWC,' )
        call output ( f_s(i_f), before=') F = ', advance='yes' )
        call dump ( p(:,:,i_theta,i_f) )
        if ( associated(dP_dIWC) .and. details >= 2 ) then
          call output ( Rad2Deg*theta_s(i_theta), before='dP_dIWC(T X IWC,' )
          call output ( f_s(i_f), before=') F = ', advance='yes' )
          call dump ( dP_dIWC(:,:,i_theta,i_f) )
          call output ( Rad2Deg*theta_s(i_theta), before='dP_dT(T X IWC,' )
          call output ( f_s(i_f), before=') F = ', advance='yes' )
          call dump ( dP_dT(:,:,i_theta,i_f) )
        end if
      end do
    end do
  end subroutine Dump_Mie

  ! ---------------------------------------------------  Read_Mie  -----
  subroutine Read_Mie ( FileName )

    use Allocate_Deallocate, only: Test_Allocate
    use IO_stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSHDF5, only: IsHDF5DSPresent, LoadPtrFromHDF5DS
    use HDF5, only: H5FClose_f, H5FOpen_f, H5F_ACC_RDONLY_F
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStrings, only: Capitalize

    character(len=*), intent(in) :: FileName

    integer :: FileID, Lun, Status

    integer :: N_f, N_IWC, N_T, N_theta, N_cut
    real(r8) :: R_min, R_max
    logical :: Derivs

    ! Destroy the old database
    call destroy_Mie

    if ( index(capitalize(fileName),'.HDF') /= 0 ) then
      ! Read tables from HDF5 file
      call H5FOpen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open Mie tables file " // trim(fileName) )
      ! Get the index arrays
      call loadPtrFromHDF5DS ( fileID, 'IWC_s', iwc_s )
      call loadPtrFromHDF5DS ( fileID, 'T_s', T_s )
      call loadPtrFromHDF5DS ( fileID, 'THETA_s', theta_s )
      call loadPtrFromHDF5DS ( fileID, 'F_s', f_s )
      ! Get the betas
      call loadPtrFromHDF5DS ( fileID, 'Beta_c_a', Beta_c_a )
      call loadPtrFromHDF5DS ( fileID, 'Beta_c_e', Beta_c_e )
      call loadPtrFromHDF5DS ( fileID, 'Beta_c_s', Beta_c_s )
      ! Get the beta derivatives
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dIWC_c_a', dBeta_dIWC_c_a )
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dT_c_a', dBeta_dT_c_a )
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dIWC_c_e', dBeta_dIWC_c_e )
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dT_c_e', dBeta_dT_c_e )
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dIWC_c_s', dBeta_dIWC_c_s )
      call loadPtrFromHDF5DS ( fileID, 'dBeta_dT_c_s', dBeta_dT_c_s )
      ! Get the phase function
      call loadPtrFromHDF5DS ( fileID, 'P', P )
      ! Get the phase function derivatives
      if ( IsHDF5DSPresent ( fileID, 'dP_dIWC' ) ) then
        call loadPtrFromHDF5DS ( fileID, 'dP_dIWC', dP_dIWC )
        call loadPtrFromHDF5DS ( fileID, 'dP_dT', dP_dT )
      end if
      call H5FClose_f ( fileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to close Mie tables file " // trim(fileName) )
    else
      ! Read tables from Fortran unformatted file
      ! Open the file
      call get_lun ( lun, msg=.false. )
      if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "No logical unit numbers available to read " // trim(fileName) )
      open ( unit=lun, file=filename, status='old', form='unformatted', &
        & access='sequential', iostat=status )
      if ( status /= 0 ) then
        call io_error ( "Unable to open Mie tables file ", status, filename )
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Unable to open Mie tables file " // trim(fileName) )
      end if

      ! Read some sizes and scalars
      read ( lun ) n_f, n_IWC, n_T, n_theta, r_min, r_max, n_cut, derivs

      ! Allocate and read the coordinate arrays
      allocate ( f_s(n_f), stat=status )
      call test_allocate ( status, moduleName, 'F_s', (/ 1 /), (/ n_f /) )
      allocate ( iwc_s(n_iwc), stat=status )
      call test_allocate ( status, moduleName, 'IWC_s', (/ 1 /), (/ n_iwc /) )
      allocate ( T_s(n_T), stat=status )
      call test_allocate ( status, moduleName, 'T_s', (/ 1 /), (/ n_T /) )
      allocate ( theta_s(n_theta), stat=status )
      call test_allocate ( status, moduleName, 'theta_s', (/ 1 /), (/ n_theta /) )
      read ( lun ) iwc_s, t_s, theta_s, f_s

      ! Allocate the phase array and its derivatives
      allocate ( p(n_t, n_iwc, n_theta, n_f), stat=status )
      call test_allocate ( status, moduleName, 'P', (/ 1,1,1,1 /), &
        & (/ n_t,n_iwc,n_theta,n_f /) )
      if ( derivs ) then
        allocate ( dP_dIWC(n_t, n_iwc, n_theta, n_f), stat=status )
        call test_allocate ( status, moduleName, 'dP_dIWC', (/ 1,1,1,1 /), &
          & (/ n_t,n_iwc,n_theta,n_f /) )
        allocate ( dP_dT(n_t, n_iwc, n_theta, n_f), stat=status )
        call test_allocate ( status, moduleName, 'dP_dT', (/ 1,1,1,1 /), &
          & (/ n_t,n_iwc,n_theta,n_f /) )
      end if

      ! Read everything
      call Read_Mie_Auto
      if ( derivs ) call Read_Mie_Derivs_Auto

      close ( lun )
    end if

  contains

    ! These are internal subroutines to avoid a bunch of explicit allocations
    subroutine Read_Mie_Auto

      integer, parameter :: I_c_e = 1, I_c_s = 2

!     real(r8) :: Beta(n_t,n_iwc,n_f,2) ! Final dimension: 1 = Beta(c_e), 2 = Beta(c_s)
!     real(r8) :: Eest(n_t,n_iwc,n_f,2) ! Error estimate
!     integer :: NFunc(2,n_t,n_iwc,n_f,2)
!     integer :: MaxOrd(n_t,n_iwc,n_f,2)
      real(r8) :: E_P(n_t,n_iwc,n_theta,n_f) ! Error estimate
      integer :: NFuncP(2, n_t, n_iwc, n_theta, n_f, 3)
      integer :: MaxOrdP(n_t, n_iwc, n_theta, n_f, 3)

!     read ( lun ) beta(:,:,:,i_c_e), eest(:,:,:,i_c_e), nFunc(:,:,:,:,i_c_e), maxOrd(:,:,:,i_c_e)
!     read ( lun ) beta(:,:,:,i_c_s), eest(:,:,:,i_c_s), nFunc(:,:,:,:,i_c_s), maxOrd(:,:,:,i_c_s)
      read ( lun ) ! Skip beta(:,:,:,i_c_e) etc
      read ( lun ) ! Skip beta(:,:,:,i_c_s) etc
      read ( lun ) p, e_p, nFuncP(:,:,:,:,:,1), maxOrdP(:,:,:,:,1)
    end subroutine Read_Mie_Auto

    subroutine Read_Mie_Derivs_Auto

      integer, parameter :: I_c_e = 1, I_c_s = 2

!     real(r8) :: dBeta_dIWC(n_t, n_iwc, n_f, 2)
!     real(r8) :: E_dBeta_dIWC(n_t, n_iwc, n_f, 2) ! Error estimate
!     real(r8) :: dBeta_dT(n_t, n_iwc, n_f, 2)
!     real(r8) :: E_dBeta_dT(n_t, n_iwc, n_f, 2) ! Error estimate
!     integer :: NFunc(2,n_t,n_iwc,n_f,3:6)
!     integer :: MaxOrd(n_t,n_iwc,n_f,3:6)
      real(r8) :: e_dP_dIWC(n_t,n_iwc,n_theta,n_f) ! Error estimate
      integer :: NFuncP(2,n_t,n_iwc,n_theta,n_f,3)
      integer :: MaxOrdP(n_t,n_iwc,n_theta,n_f,3)
      real(r8) :: e_dP_dT(n_t,n_iwc,n_theta,n_f) ! Error estimate

!     read ( lun ) dBeta_dIWC(:,:,:,i_c_e),  e_dBeta_dIWC(:,:,:,i_c_e), nFunc(:,:,:,:,3), maxOrd(:,:,:,3)
!     read ( lun ) dBeta_dIWC(:,:,:,i_c_s),  e_dBeta_dIWC(:,:,:,i_c_s), nFunc(:,:,:,:,4), maxOrd(:,:,:,4)
!     read ( lun ) dBeta_dT(:,:,:,i_c_e),  e_dBeta_dT(:,:,:,i_c_e), nFunc(:,:,:,:,5), maxOrd(:,:,:,5)
!     read ( lun ) dBeta_dT(:,:,:,i_c_s),  e_dBeta_dT(:,:,:,i_c_s), nFunc(:,:,:,:,6), maxOrd(:,:,:,6)
      read ( lun ) ! skip dBeta_dIWC(:,:,:,i_c_e) etc
      read ( lun ) ! skip dBeta_dIWC(:,:,:,i_c_s) etc
      read ( lun ) ! skip dBeta_dT(:,:,:,i_c_e) etc
      read ( lun ) ! skip dBeta_dT(:,:,:,i_c_s) etc
      read ( lun ) dP_dIWC, e_dP_dIWC, nFuncP(:,:,:,:,:,2), maxOrdP(:,:,:,:,2)
      read ( lun ) dP_dT, e_dP_dT, nFuncP(:,:,:,:,:,3), maxOrdP(:,:,:,:,3)
    end subroutine Read_Mie_Derivs_Auto

  end subroutine Read_Mie

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Read_Mie_m

! $Log$
! Revision 2.2  2008/06/05 02:15:26  vsnyder
! Added HDF, spiffed up dump
!
! Revision 2.1  2008/05/20 00:26:40  vsnyder
! Initial commit
!
