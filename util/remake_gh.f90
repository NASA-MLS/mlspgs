program REMAKE_GH

  implicit NONE

  real, parameter :: DTEMOD(13) = (/ 1945., 1950., 1955., 1960., 1965., &
    &                                1970., 1975., 1980., 1985., 1990., &
    &                                1995., 2000., 2005. /)

  real :: ERAD

  character*12, parameter :: FILMOD(13) = (/ 'dgrf45.dat  ', 'dgrf50.dat  ', &
    &                        'dgrf55.dat  ', 'dgrf60.dat  ', 'dgrf65.dat  ', &
    &                        'dgrf70.dat  ', 'dgrf75.dat  ', 'dgrf80.dat  ', &
    &                        'dgrf85.dat  ', 'dgrf90.dat  ', 'dgrf95.dat  ', &
    &                        'igrf00.dat  ', 'igrf00s.dat ' /)

  real :: GH(144)

  integer :: I, IER, J, K
  integer, parameter :: IU = 10
  integer :: NGH, NGHMAX = 0, NMAX, NMAXMAX = 0
  character(len=127) :: PREFIX

  call getarg ( 1, prefix )
  if ( prefix == '' ) then
    print *, 'Enter prefix of file names: '
    read ( *, '(a)' ) prefix
  end if

  open ( 99, status='scratch', form='unformatted' )

  do i = 1, size(filmod)

    call getshc ( iu, trim(prefix) // filmod(i), nmax, gh, ier, erad, ngh )

    write ( 99 ) nmax, ngh, erad, gh(:ngh)

    nmaxmax = max(nmaxmax, nmax)
    nghmax = max(nghmax, ngh)
  end do
  rewind ( 99 )

  write ( *, '(a)' ) '  type :: GH_T', &
    &                '    real    :: YEAR', &
    &                '    character(len=12) :: FILE', &
    &                '    integer :: NMAX', &
    &                '    integer :: NGH', &
    &                '    real ::    ERAD'
  write ( *, '(a:,i3,a)' ) &
    &                '    real ::    GH(', nghmax, ')', &
    &                '  end type GH_T'
  write ( *, '(a,i2,a)' ) &
    &                '  type(gh_t), save :: GHS(', size(filmod), ')'

  do i = 1, size(filmod)
    read ( 99 ) nmax, ngh, erad, gh(:ngh)
    write ( *, '(a,i3,a)' ) '  data GHS(', i, ') / gh_t ( &'
    write ( *, '(a,f7.1,3a,2(i3,a),g14.7,a)' ) &
      & '    & ', dtemod(i), ', "', filmod(i), '", ', &
      &           nmax, ', ', ngh, ', ', erad, ', (/ &'
    gh(ngh+1:) = 0.0
    do j = 1, nghmax, 5
      k = min(j+4,nghmax)
      write ( *, '("    & ",5(g14.7:,", "))', advance='no' ) gh(j:k)
      if ( k /= nghmax ) then
        write ( *, '(", &")' )
      else
        write ( *, '("/) ) /")' )
      end if
    end do
  end do
  close ( 99 )

contains

  ! -----------------------------------------------------  GETSHC  -----
  subroutine GETSHC ( IU, FSPEC, NMAX, GH, IER, ERAD, NGH )

! ===============================================================
!
!       Version 1.01
!
!       Reads spherical harmonic coefficients from the specified
!       file into an array.
!
!       Input:
!           IU     - Logical unit number
!           FSPEC  - File specification
!
!       Output:
!           NMAX   - Maximum degree and order of model
!           GH     - Schmidt quasi-normal internal spherical
!                    harmonic coefficients
!           IER    - Error number: =  0, no error
!                                  = -2, records out of order
!                                  = FORTRAN run-time error number
!           ERAD   - Earth's radius associated with the spherical
!                    harmonic coefficients, in the same units as
!                    elevation.
!           NGH    - Number of elements in GH
!
!       A. Zunde
!       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================

    integer, intent(in) :: IU
    character(len=*), intent(in) :: FSPEC
    integer, intent(out) :: NMAX
    real, intent(out) :: GH(:)
    integer, intent(out) :: IER
    real, intent(out) :: ERAD
    integer, intent(out), optional :: NGH

    real :: G, H
    integer :: I, N, NN, M, MM

! ---------------------------------------------------------------
!       Open coefficient file. Read past first header record.
!       Read degree and order of model and Earth's radius.
! ---------------------------------------------------------------
    open ( iu, file=fspec, status='old', iostat=ier, err=999 )

    read ( iu, *, iostat=ier, err=999 ) ! Skip the file name in the file
    read ( iu, *, iostat=ier, err=999 ) nmax, erad
! ---------------------------------------------------------------
!       Read the coefficient file, arranged as follows:
!
!                                       N     M     G     H
!                                       ----------------------
!                                   /   1     0    GH(1)  -
!                                  /    1     1    GH(2) GH(3)
!                                 /     2     0    GH(4)  -
!                                /      2     1    GH(5) GH(6)
!           NMAX*(NMAX+3)/2     /       2     2    GH(7) GH(8)
!              records          \       3     0    GH(9)  -
!                                \      .     .     .     .
!                                 \     .     .     .     .
!           NMAX*(NMAX+2)          \    .     .     .     .
!           elements in GH          \  NMAX  NMAX   .     .
!
!       N and M are, respectively, the degree and order of the
!       coefficient.
! ---------------------------------------------------------------

    i = 0                                                   
o:  do nn = 1, nmax                                         
      do mm = 0, nn                                       
        read (iu, *, iostat=ier, err=999, end=999) n, m, g, h
        if ( nn /= n .or. mm /= m ) then                
          ier = -2                                    
          exit o
        end if                                          
        i = i + 1                                       
        gh(i) = g                                       
        if ( m /= 0 ) then                              
          i = i + 1                                   
          gh(i) = h                                   
        end if                                          
      end do                                                 
    end do o

999 close (iu)
    if ( present(ngh) ) ngh = i

  end subroutine GETSHC

end program REMAKE_GH
