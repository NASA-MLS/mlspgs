! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! ifort-v17 -o Test_PGS_CSC_TeoToECR -Bstatic Test_PGS_CSC_GeoToECR.f90
! -I ../../lib/IFC.Linux_17 -I /software/toolkit/ifc17/toolkit/include
! -L ../../lib/IFC.Linux_17 -lmls -L ../../bin/IFC.Linux_17 -lutctotai
! -L/software/toolkit/ifc17/toolkit/lib/linux -lPGSTK
! -L/software/toolkit/ifc17/hdfeos/lib/linux -lhdfeos
! -L/software/toolkit/ifc17/hdfeos5/lib/linux -lhe5_hdfeos
! -L/software/toolkit/ifc17/hdf/lib -lmfhdf -ldf -lz -ljpeg -lsz
! -L/software/toolkit/ifc17/hdfeos/lib/linux -lGctp

program ECRtoECI

  ! Read the Earth ellipse (default WGS84), the time, and ECR
  ! coordinates, from the command line.  Print ECI coordinates.
  ! Optionally, accept a -I argument to specify a file from which to
  ! read ECR coordinates (time and ellipse still on the command line).

  use, intrinsic :: ISO_Fortran_ENV, only: Input_Unit
  use SDPToolkit, only: PGS_CSC_ECRtoECI, PGS_S_Success, PGS_SMF_GetMsg

  implicit NONE

  character(28) :: AsciiUTC
  character (25) :: EarthEllipseTag = 'WGS84'
  double precision :: ECI(6) = 0, ECR(6) = 0, Offsets(1) = 0.0d0
  logical :: File = .false.
  integer :: I, J
  character(127) :: IOMSG
  character(127) :: Line ! Command-line argument, or from the file
  integer :: LUN = 10    ! For input
  integer :: Stat

  call get_environment_variable ( 'PGSHOME', line )
  if ( line == '' ) then
    print '(a)', 'The environment variable PGSHOME is not set'
    stop
  end if

  i = 0
o:do
    i = i + 1
    call get_command_argument ( i, line )
    select case ( line(1:2) )
    case ( '-e' )
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      earthEllipseTag = adjustl(line(3:))
    case ( '-h' )
      call get_command_argument ( 0, line )
      print '(3a)', 'Usage: ', trim(line), ' [options] Time [3 ECR coordinates]'
      print '(a)',  ' Time is yyyy-dddThh:mm:ss.sss'
      print '(a)',  ' Options:'
      print '(2a)', '  -e[ ]tag  => Earth Ellipse Tag, default ', trim(earthEllipseTag)
      print '(a)',  '  -h        => This output'
      print '(a)',  '  -i[ ]file => Read more input from file, -i - for stdin'
      stop
    case ( '-I', '-i' )
      if ( line(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, line(3:) )
      end if
      file = .true.
      if ( adjustl(line(3:)) == '-' ) then
        lun = input_unit
      else
        open ( 10, file=adjustl(line(3:)), status='old', iostat=stat, iomsg=iomsg )
        if ( stat /= 0 ) then
          print '(a,i0,2a/a)', 'Error ', stat, ' while trying to open ', &
            & trim(adjustl(line(3:))), trim(iomsg)
          stop
        end if
      end if
    case default
      asciiUTC = line
      do j = 1, 3
        i = i + 1
        call get_command_argument ( i, line )
        if ( line == '' ) then
          if ( file .and. j == 1 ) exit o
          print '(a)', 'Too few coordinates given'
          stop
        end if
        read ( line, *, iostat=stat, iomsg=iomsg ) ecr(j)
        if ( stat /= 0 ) then
          print '(a,i0,2a/a)', 'Error ', stat, ' while trying to read ', &
            & trim(adjustl(line(3:))), trim(iomsg)
          stop
        end if
      end do
      exit
    end select
  end do o

  do
    if ( line /= '' ) then
      line = asciiUTC
      stat = PGS_CSC_ECRtoECI (1, line, offsets, ECR, ECI)
      if ( stat /= PGS_S_Success ) then
        print '(a,i0)', 'Status from PGS_CSC_ECRtoECI = ', stat
        print '(a)', 'Is the environment variable PGSHOME set?'
        call PGS_SMF_GetMsg ( stat, line, iomsg )
        print '(a)', trim(line), trim(iomsg)
        stop
      end if
      print '(1p,3g15.7)', ECI(1:3)
    end if
    if ( .not. file ) exit
    read ( lun, '(a)', end=9 ) line
    read ( line, *, iostat=stat, iomsg=iomsg ) ECR(1:3)
    if ( stat /= 0 ) then
      print '(a,i0,2a/a)', 'Error ', stat, ' while trying to read ', &
        & trim(adjustl(line(3:))), trim(iomsg)
      stop
    end if
  end do

9 continue

contains

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = IdParm == ''
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end program ECRtoECI

! $Log$
! Revision 1.1  2017/01/05 00:29:15  vsnyder
! Initial commit
!
