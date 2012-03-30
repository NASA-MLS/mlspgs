! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program REMAKE_GH

! Read spherical harmonic coefficients for Earth's magnetic field models
! and write a type definition, a variable declaration, and data statements
! for them, as a Fortran include file.

! Command line arguments:
!   -p prefix  => a prefix for input file names, default "."
!   -o output  => output file, default stdout
!   files...   => up to 20 input file names

  use, intrinsic :: ISO_Fortran_Env, only: Output_Unit

  implicit NONE

  character(127) :: ARG
  integer :: DT(8)  ! Output from date_and_time intrinsic subroutine
  real :: DTEMOD(20)
  real :: ERAD
  character(13) :: FILMOD(20) = '*************'
  real :: GH(1023)
  integer :: I, IER, J, K
  integer :: NFILE, NGH, NGHMAX = 0, NMAX, NMAXMAX = 0, OU = output_unit
  character(len=127) :: OFILE = '', PREFIX = ''

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  i = 0
  do
    i = i + 1
    call get_command_argument ( i, arg )
    select case ( arg(1:3) )
    case ( '-o ' )
      i = i + 1
      call get_command_argument ( i, ofile )
      ou = 42
      open ( ou, file=trim(ofile), form='formatted' )
    case ( '-p ' )
      i = i + 1
      call get_command_argument ( i, prefix )
    case ( '-h ' )
      call get_command_argument ( 0, arg )
      print '(3a)', 'Usage: ', trim(arg), ' [options] input-file-names'
      print '(a)' , ' Options: -o file => output file name, else stdout'
      print '(a)' , '          -p prefix => prefix for input file names'
      print '(a)' , '                       / is appended to prefix'
      stop
    case default
      exit
    end select
  end do

  nfile = 0
  do while ( i <= command_argument_count() )
    nfile = nfile + 1
    call get_command_argument ( i, filmod(nfile) )
    i = i + 1
  end do

  if ( prefix == '' ) prefix = '.'

  call date_and_time ( values=dt )
  write ( ou, '(a,i4.4,2("-",i2.2),"T",i2.2,2(":",i2.2),".",i0,2a)' ) &
    & '! Created at ', dt(1:3), dt(5:8), ' from files in ', trim(prefix)

  open ( 99, status='scratch', form='unformatted' )

  ! Read the files
  do i = 1, nfile

    call getshc ( trim(prefix) // '/' // trim(filmod(i)), nmax, dtemod(i), gh, ier, erad, ngh )

    write ( 99 ) nmax, ngh, erad, gh(:ngh)

    nmaxmax = max(nmaxmax, nmax)
    nghmax = max(nghmax, ngh)
  end do
  rewind ( 99 )

  ! Write a type definition for the data
  write ( ou, '(a)' ) '  type :: GH_T', &
    &                 '    real    :: YEAR', &
    &                 '    character(len=12) :: FILE', &
    &                 '    integer :: NMAX', &
    &                 '    integer :: NGH', &
    &                 '    real ::    ERAD'
  write ( ou, '(a:,i0,a)' ) &
    &                 '    real ::    GH(', nghmax, ')', &
    &                 '  end type GH_T'
  ! Write a declaration for the variable
  write ( ou, '(a,i0,a)' ) &
    &                 '  type(gh_t), save :: GHS(', nfile, ')'

  ! Write the data statements
  do i = 1, nfile
    read ( 99 ) nmax, ngh, erad, gh(:ngh)
    write ( ou, '(a,i0,a)' ) '  data GHS(', i, ') / gh_t ( &'
    write ( ou, '(a,f7.1,3a,2(i3,a),g14.7,a)' ) &
      & '    & ', dtemod(i), ', "', filmod(i), '", ', &
      &           nmax, ', ', ngh, ', ', erad, ', (/ &'
    gh(ngh+1:) = 0.0
    do j = 1, nghmax, 5
      k = min(j+4,nghmax)
      write ( ou, '("    & ",5(g14.7:,", "))', advance='no' ) gh(j:k)
      if ( k /= nghmax ) then
        write ( ou, '(", &")' )
      else
        write ( ou, '("/) ) /")' )
      end if
    end do
  end do
  close ( 99 )

contains

  ! -----------------------------------------------------  GETSHC  -----
  subroutine GETSHC ( FSPEC, NMAX, YEAR, GH, IER, ERAD, NGH )

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

    character(len=*), intent(in) :: FSPEC
    integer, intent(out) :: NMAX
    real, intent(out) :: YEAR
    real, intent(out) :: GH(:)
    integer, intent(out) :: IER
    real, intent(out) :: ERAD
    integer, intent(out), optional :: NGH

    integer, parameter :: IU = 10

    real :: G, H
    integer :: I, N, NN, M, MM

! ---------------------------------------------------------------
!       Open coefficient file. Read past first header record.
!       Read degree and order of model and Earth's radius.
! ---------------------------------------------------------------

    open ( iu, file=fspec, status='old', iostat=ier, err=666 )

    read ( iu, *, iostat=ier, err=666 ) ! Skip the file name in the file
    read ( iu, *, iostat=ier, err=666 ) nmax, erad, year

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
        read (iu, *, iostat=ier, err=666, end=999) n, m, g, h
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

    return

666 print '(a,i0/2a)', 'I/O Error, iostat = ', ier, 'File = ', trim(fspec)
  end subroutine GETSHC

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
!---------------------------------------------------------------------------

end program REMAKE_GH

! $Log$
! Revision 1.7  2012/03/29 20:21:24  vsnyder
! Add -h option for usage summary
!
! Revision 1.6  2012/03/29 00:43:08  vsnyder
! Add time stamp, make GH array bigger
!
! Revision 1.5  2012/03/14 00:58:03  vsnyder
! Remove goofy method to determine number of files
!
! Revision 1.4  2012/03/14 00:57:08  vsnyder
! Add copyright statement
!
! Revision 1.3  2012/03/14 00:54:01  vsnyder
! Revised to get file names from command line
!
! Revision 1.2  2012/03/13 23:46:57  vsnyder
! Revise to get file names from input instead of internal data
!
! Revision 1.1  2012/03/13 23:15:55  vsnyder
! Initial commit
!
