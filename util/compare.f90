! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program COMPARE

! Compare two files and report the maximum relative and absolute differences
! between the numbers in them.

! In each file, look for a line that ends with \ <number>.

! If the lines are different, report an error and stop.

! Read <number> numbers from each of the two files.  Compute their maximum
! relative and absolute difference. If this is different from zero and the
! -a option is specified, print the line with the \ and the difference.

! Repeat to the end of the files.

! Print the maximum relative and absolute difference anywhere at the end
! if it's not zero.

  use Machine, only: IO_Error

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  !---------------------------------------------------------------------------

  real, allocatable, dimension(:) :: AD      ! Absolute difference of R1, R2
  logical :: All = .false.   ! Show nonzero differences for all quantities
  real :: AMAX               ! Maximum absolute difference for one R1, R2 pair
  real :: AMAXG = 0.0        ! Global maximum of all values of AMAX
  logical :: AnyNaN(2) = .false. ! R1 or R2 has a NaN
  logical :: CONT = .false.  ! Continue even if control lines differ
  logical :: END
  character(127) :: File1, File2
  integer :: I, J
  integer :: LAMAX, LRMAX    ! Locations of absolute, relative max. diffs.
  logical :: LOUD = .true.   ! Messages about unequal file lengths etc.
  character(127) :: Line1, line2
  integer :: N
  real, allocatable, dimension(:) :: R1, R2  ! Inputs
  real, allocatable, dimension(:) :: RD      ! Relative difference of R1, R2
  real :: RMAX               ! Maximum relative difference for one R1, R2 pair
  real :: RMAXE = -huge(0.0) ! RMAX in units of epsilon
  real :: RMAXG = -huge(0.0) ! Global maximum of all values of RMAX
  integer :: Status
  logical :: Zero = .false.  ! If ( all ), show zero differences, too.

  i = 1
  do
    call getarg ( i, line1 )
    if ( line1(1:1) /= '-' .and. line1 /= '' ) exit
    if ( line1(1:2) == '- ' ) then
      i = i + 1
      exit
    end if
    if ( line1(1:1) == '-' ) then
      do j = 2, len_trim(line1)
        if ( line1(j:j) == 'a' ) then
          all = .true.
        else if ( line1(j:j) == 'c' ) then
          cont = .true.
        else if ( line1(j:j) == 'q' ) then
          loud = .false.
        else if ( line1(j:j) == 'z' ) then
          zero = .true.
        else
          call usage
        end if
      end do
    else
      call usage
    end if
    i = i + 1
  end do

  call getarg ( i, file1 )
  open ( 10, file=file1, form='formatted', status='old', iostat=status )
  if ( status /= 0 ) then
    call io_error ( 'Unable to open input file', status, file1 )
    stop
  end if
  call getarg ( i+1, file2 )
  open ( 11, file=file2, form='formatted', status='old', iostat=status )
  if ( status /= 0 ) then
    call io_error ( 'Unable to open input file', status, file2 )
    stop
  end if

  if ( all ) then
    print *, '  Absolute           Relative            Relative'
    print *, '  Difference     At  Difference      At  Epsilons      After'
  end if

  do

    do
      read ( 10, '(a)', iostat=status ) line1
      end = status /= 0
      if ( end ) exit
      if ( index(line1,'\') /= 0 ) exit
    end do

    do
      read ( 11, '(a)', iostat=status ) line2
      if ( status /= 0 ) exit
      i = index(line2,'\')
      if ( i /= 0 ) exit
    end do

    if ( (end .neqv. status /= 0) .and. loud )  print *, 'Input file lengths unequal'
    if ( end .or. status /= 0 ) exit

    if ( line1 /= line2 ) then
      if ( loud ) then
        print *, 'Control lines differ:'
        print *, trim(line1)
        print *, trim(line2)
      end if
      if ( .not. cont ) exit
    end if

    read ( line2(i+1:), *, iostat=status ) n
    if ( status /= 0 ) then
      call io_error ( line2, status )
      exit
    end if

    allocate ( ad(n), r1(n), r2(n), rd(n) )

    read ( 10, *, iostat=status ) r1
    if ( status /= 0 ) exit
    read ( 11, *, iostat=status ) r2
    if ( status /= 0 ) exit

    ! DO NOT apply deMauvre's theorem.  If you write
    !    any ( r1 > 0.0 .and. r1 < 0.0 ) it won't catch NaN's
    if ( any( .not. (r1 <= 0.0 .or. r1 >= 0.0) ) ) anyNaN(1) = .true.
    if ( any( .not. (r2 <= 0.0 .or. r2 >= 0.0) ) ) anyNaN(2) = .true.

    ad = abs(r1 - r2)
    lamax = maxloc(ad,1)
    amax = ad(lamax)
    if ( amax > 0.0 .or. all .and. zero ) then
      rd = 2.0 * ad / abs(r1 + r2)
      lrmax = maxloc(rd,1)
      rmax = rd(lrmax)
      rmaxe = rmax / epsilon(rd)
      if ( all ) then
        print '(1pg14.8,i6,1pg14.8,i6,1pg15.8,1x,a)', &
          & amax, lamax, rmax, lrmax, rmaxe, trim(line1)
!       print *, 'After ', trim(line1), ', Maximum difference =', amax, ' at', lamax
!       print *, 'Relative =', rmax, ' =', rmaxe, ' epsilons', ' at', lrmax
      end if
    end if
    rmaxg = max(rmaxe,rmaxg)
    amaxg = max(amax,amaxg)

    deallocate ( ad, r1, r2, rd )

  end do

  if ( rmaxg > 0.0 .or. zero ) then
    print *, 'Maximum relative difference anywhere =', rmaxg, 'epsilons =', &
      & rmaxg*epsilon(rmaxg)
    print *, 'Maximum absolute difference anywhere =', amaxg
  end if

  if ( anyNaN(1) ) print *, trim(file1), ' has a NaN somewhere'
  if ( anyNaN(2) ) print *, trim(file2), ' has a NaN somewhere'

contains
  subroutine USAGE
    call getarg ( 0, line1 )
    print *, 'Usage: ', trim(line1), ' [option] file1 file2'
    print *, ' Options: -a => Show nonzero difference for all quantities'
    print *, '          -c => Continue even if control lines differ'
    print *, '          -q => No messages about unequal file lengths etc.'
    print *, '          -z => Show zero difference summary at the end'
    print *, '                With -a, show zero individual differences too.'
    stop
  end subroutine USAGE

end program

! $Log$
! Revision 1.1  2003/07/03 18:08:34  vsnyder
! Initial commit
!
