! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

  implicit NONE

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  integer, parameter :: RK = kind(0.0d0) ! kind for calculations
  integer, parameter :: RS = kind(0.0e0) ! kind for * output
  real(rk), allocatable, dimension(:) :: AD      ! Absolute difference of R1, R2
  logical :: All = .false.        ! Show nonzero differences for all quantities
  real(rk) :: AMAX                ! Maximum absolute difference for one R1, R2 pair
  real(rk) :: AMAXG = 0.0         ! Global maximum of all values of AMAX
  real(rk) :: AbsAtRmaxG          ! Absolute error at maximum relative error
  logical :: AnyNaN(3) = .false.  ! R1, R2 or a relative difference has a NaN
  real(rk) :: AVG(2), AVGSR(2), AVGSA(2)
  logical :: CONT = .false.       ! Continue even if control lines differ
  logical :: DoStats = .false.    ! -s option specified
  logical :: END
  character(127) :: Errmsg
  logical :: Error
  character(127) :: File1, File2
  integer :: I, I1, I2, J
  logical :: Same = .false.       ! Print "Identical" if files are the same
  integer :: LAMAX, LRMAX         ! Locations of absolute, relative max. diffs.
  logical :: LOUD = .true.        ! Messages about unequal file lengths etc.
  character(127) :: Line1, line2
  integer :: Loc1, Loc2           ! Line numbers in files
  logical, allocatable, dimension(:) :: M    ! abs(r1+r2) > 0.0
  integer :: N, N1, N2            ! Numbers of values
  character(len=3) :: NaNString
  real(rk) :: NaNValue
  logical :: NoZero = .false.     ! Don't compare where either one is zero, -Z
  real(rk), allocatable, dimension(:) :: R1, R2  ! Inputs
  real(rk), allocatable, dimension(:) :: RD      ! Relative difference of R1, R2
  real(rk) :: RelAtAmaxG
  real(rk) :: RMAX = -huge(0.0)   ! Maximum relative difference for one R1, R2 pair
  character(127) :: RMAXB = ''    ! Label of block having largest value of RMAXV
  real(rk) :: RMAXG = -huge(0.0)  ! Global maximum of all values of RMAX
  real(rk) :: RMAXL = -huge(0.0)  ! Relative at max abs diff (LAMAX)
  real(rk) :: RMAXV = -huge(0.0)  ! Difference relative to VMAX
  real(rk) :: RMAXVG = -huge(0.0) ! Global maximum of all values of RMAXV
  integer :: Status
  real(rk) :: STDEV(2), STDEVR(2), STDEVA(2)
  logical :: Verbose = .false.
  real(rk) :: VMAX                ! Maximum absolute value in R1 or R2
  logical :: Zero = .false.       ! If ( all ), show zero differences, too.

  NaNString = "NaN"
  read ( NaNString, * ) NaNvalue
  ! read ( (/("NaN", i = 1, 10 ) /), * ) AbsAtRmaxG, RelAtAmaxG, avgsr, avgsa, &
  !  & stdevr, stdeva
  AbsAtRmaxG = NaNvalue
  RelAtAmaxG = NaNvalue
  avgsr = NaNvalue
  avgsa = NaNvalue
  same = .true. ! in case files are identical, even without -i option
  stdevr = NaNvalue
  stdeva = NaNvalue

  i = 1
  do
    call get_command_argument ( i, line1 )
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
        else if ( line1(j:j) == 'i' ) then
          same = .true.
        else if ( line1(j:j) == 'q' ) then
          loud = .false.
        else if ( line1(j:j) == 's' ) then
          doStats = .true.
        else if ( line1(j:j) == 'v' ) then
          print *, Id
        else if ( line1(j:j) == 'V' ) then
          verbose = .true.
        else if ( line1(j:j) == 'z' ) then
          zero = .true.
        else if ( line1(j:j) == 'Z' ) then
          noZero = .true.
        else
          call usage
        end if
      end do
    else
      call usage
    end if
    i = i + 1
  end do

  error = .false.
  call get_command_argument ( i, file1 )
  open ( 10, file=file1, form='formatted', status='old', iostat=status, &
    & iomsg=errmsg )
  if ( status /= 0 ) then
    print 1, trim(file1), i, status, trim(errmsg)
1   format ( 'Unable to open input file "', a, '" in argument ', i0 / &
      & 'Status = ', i0, ', Message = ', a )
    error = .true.
  end if
  call get_command_argument ( i+1, file2 )
  open ( 11, file=file2, form='formatted', status='old', iostat=status, &
    & iomsg=errmsg )
  if ( status /= 0 ) then
    print 1, trim(file2), i+1, status, trim(errmsg)
    error = .true.
  end if
  if ( error ) stop

  if ( all ) then
    print '(a)', '   Max Abs     Max Abs            Rel at Max  Max Rel            Abs at Max  Rel to Max'
    print '(a)', '   Value       Diff          At   Abs Diff    Diff          At   Rel Diff    Abs Value'
  end if

  loc1 = 0
  loc2 = 0

  do

    do
      loc1 = loc1 + 1
      read ( 10, '(a)', iostat=status ) line1
      end = status /= 0
      if ( end ) exit
      i1 = index(line1,'\')
      if ( i1 /= 0 ) exit
    end do

    do
      loc2 = loc2 + 1
      read ( 11, '(a)', iostat=status ) line2
      if ( status /= 0 ) exit
      i2 = index(line2,'\')
      if ( i2 /= 0 ) exit
    end do

    if ( verbose ) print *, 'Read control lines from both files'
    if ( (end .neqv. status /= 0) ) then
      if ( loud ) print '(a)', 'Input file lengths unequal'
      same = .false.
    end if
    if ( end .or. status /= 0 ) exit

    if ( line1 /= line2 ) then
      if ( loud ) then
        print '(a)', 'Control lines differ:'
        print '(i5,2a)', loc1, ': ', trim(line1)
        print '(i5,2a)', loc2, ': ', trim(line2)
      end if
      if ( .not. cont ) exit
    end if

    read ( line1(i1+1:), *, iostat=status, iomsg=errmsg ) n1
    if ( status /= 0 ) then
      print '(a,i0,2a)', 'Unable to read number N from ', loc1, ': ', &
        & trim(line1(i1+1:))
      print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
      exit
    end if
    read ( line2(i2+1:), *, iostat=status, iomsg=errmsg ) n2
    if ( status /= 0 ) then
      print '(a,i0,2a)', 'Unable to read number N from ', loc2, ': ', &
        & trim(line2(i2+1:))
      print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
      exit
    end if

    if ( n1 /= n2 ) then
      print '(2(a,i0),a,2(/i5,2a))', 'Block sizes ', n1, ' and ', n2, ' differ', &
        & loc1, ': ', trim(line1(i1+1:)), loc2, ': ', trim(line2(i2+1:))
      exit
    end if

    n = n1
    if ( verbose ) print *, 'Allocating space for blocks of ', n
    allocate ( ad(n), r1(n), r2(n), rd(n), m(n) )

    if ( n == 1 ) then
      read ( line1(i1+1:), *, iostat=status ) n, r1
      if ( status /= 0 ) then
        print '(a,i0,2a)', 'Unable to read numbers N, R1 from ', loc1, &
          & ': ', trim(line1(i1+1:))
        print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
        exit
      end if
      read ( line2(i2+1:), *, iostat=status ) n, r2
      if ( status /= 0 ) then
        print '(a,i0,2a)', 'Unable to read numbers N, R1 from ', loc2, &
          & ': ', trim(line2(i2+1:))
        print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
        exit
      end if
    else
      loc1 = loc1 + 1
      read ( 10, *, iostat=status ) r1
      if ( status /= 0 ) then
        print '(a,i0)', 'Unable to read number R1 from first input file at line ', &
          & loc1
        print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
        exit
      end if
      loc2 = loc2 + 1
      read ( 11, *, iostat=status ) r2
      if ( status /= 0 ) then
        print '(a,i0)', 'Unable to read number R1 from second input file at line ', &
          & loc2
        print '("Status = ", i0, ", Message = ",a)', status, trim(errmsg)
        exit
      end if
    end if

    if ( verbose ) print *, 'Read blocks from both files of size', size(r1)

    ! DO NOT apply deMauvre's theorem.  If you write
    !    any ( r1 > 0.0 .and. r1 < 0.0 ) it won't catch NaN's
    if ( any( .not. (r1 <= 0.0 .or. r1 >= 0.0) ) ) anyNaN(1) = .true.
    if ( any( .not. (r2 <= 0.0 .or. r2 >= 0.0) ) ) anyNaN(2) = .true.

    ad = abs(r1 - r2)
    if ( noZero ) then
      where ( r1 == 0 .or. r2 == 0 ) ad = 0
    end if
    lamax = maxloc(ad,1)
    amax = ad(lamax)
    if ( amax > 0.0 .or. all .and. zero .or. doStats ) then
      vmax = max(maxval(abs(r1)), maxval(abs(r2)))
      rd = abs(r1 + r2)
      m = rd > 0.0
      rd = 2.0 * ad / rd
      if ( any(m) ) then
        lrmax = maxloc(rd,1,mask=m)
        rmax = rd(lrmax)
      else
        lrmax = -1
        rmax = 0.0
      end if
      if ( vmax > 0.0 ) then
        rmaxv = amax / vmax
      else if ( amax > 0.0 ) then
        rmaxv = -huge(0.0)
      else
        rmaxv = 0.0
      end if
      if ( .not. ( rmax <= 0.0 .or. rmax >= 0.0 ) ) anyNaN(3) = .true.
      if ( all ) then
        rmaxl = 0
        if ( abs(r1(lamax)+r2(lamax)) > 0 ) &
          & rmaxl = 2.0 * abs(r1(lamax)-r2(lamax)) / abs(r1(lamax)+r2(lamax))
        print '(1p,2(2g12.5,i7),2g12.5,1x,a)', vmax, &
          & amax, lamax, rmaxl, &
          & rmax, lrmax, 2.0 * abs(r1(lrmax)-r2(lrmax)), rmaxv, trim(line1)
      end if
      if ( doStats ) then
        call stats ( r1, avg(1), stdev(1) )
        call stats ( r2, avg(2), stdev(2) )
        if ( all ) print '(a,1p,2g14.7,a,2g14.7)', 'Averages =', avg, ' Std. Devs. =', stdev
      end if
    end if
    if ( rmax > rmaxg ) then
      rmaxg = rmax
      absAtRmaxG = amax
      if ( doStats ) then
        avgsr = avg
        stdevr = stdev
      end if
    end if
    if ( amax > amaxg ) then
      amaxg = amax
      relAtAmaxG = rmax
      if ( doStats ) then
        avgsa = avg
        stdeva = stdev
      end if
    end if
    if ( rmaxv > rmaxvg ) then
      rmaxvg = rmaxv
      rmaxb = line1
    end if

    deallocate ( ad, r1, r2, rd, m )

  end do

  if ( rmaxvg > 0.0 .or. zero .or. anyNan(3) ) then
    print '(a/1p,5g13.5,1x,a)', &
      & " RelMaxVal    MaxRel       where MaxAbs MaxAbs       where MaxRel block", &
      & rmaxvg, rmaxg, absAtRmaxG, amaxG, relAtAmaxG, trim(rmaxb)
    if ( doStats ) &
      & print '(a/1p,8g13.6)', &
      & " Avg Rel                   Std Dev Rel               Avg Abs                   Std Dev Abs", &
      & avgsr, stdevr, avgsa, stdeva
  else if ( same ) then
    print *, 'Identical'
  else ! Print this even -q option, but there's no input at all.
    print *, 'File lengths unequal, probably one was empty'
  end if

  if ( anyNaN(1) ) print *, trim(file1), ' has a NaN somewhere'
  if ( anyNaN(2) ) print *, trim(file2), ' has a NaN somewhere'

contains

  subroutine Stats ( A, Avg, Stdev )
  ! Compute the average and standard deviation of the nonzero elements of A.
    real(rk), intent(in) :: A(:)
    real(rk), intent(out) :: Avg, Stdev
    integer :: I, N
    n = 0
    avg = 0.0
    stdev = 0.0
    do i = 1, size(A)
      if ( a(i) /= 0.0 ) then
        avg = avg + a(i)
        n = n + 1
      end if
    end do
    if ( n > 1 ) then
      avg = avg / n
      do i = 1, size(A)
        if ( a(i) /= 0.0 ) stdev = stdev + (a(i)-avg)**2
      end do
      stdev = sqrt(stdev/(n-1))
    else
      stdev = huge(0.0)
    end if
  end subroutine Stats

  subroutine USAGE
    call get_command_argument ( 0, line1 )
    print *, 'Usage: ', trim(line1), ' [option] file1 file2'
    print *, ' Options: -a => Show nonzero difference for all quantities'
    print *, '          -c => Continue even if control lines differ'
    print *, '          -i => Print "Identical" if files are identical'
    print *, '          -q => No messages about unequal file lengths etc.'
    print *, '          -s => Compute average and standard deviation of nonzero elements'
    print *, '          -v => Print the version'
    print *, '          -V => Be verbose'
    print *, '          -z => Show zero difference summary at the end'
    print *, '                With -a, show zero individual differences too.'
    print *, '          -"anything else", or missing file1 or file2'
    print *, '             => This explanation.'
    stop
  end subroutine USAGE

end program

! $Log$
! Revision 1.32  2018/05/24 03:25:39  vsnyder
! Use the same precision in all output
!
! Revision 1.31  2018/04/25 20:53:24  vsnyder
! Repair header for summary-only case
!
! Revision 1.30  2018/04/10 22:19:32  vsnyder
! Check both input files' availability before stopping if one is not
! available.  Report the argument number(s) of missing input files.
!
! Revision 1.29  2018/04/05 01:34:34  vsnyder
! Don't claim file sizes are different if files are identical
!
! Revision 1.28  2017/11/29 01:35:44  vsnyder
! Don't print 'Identical' if a file was empty
!
! Revision 1.27  2017/11/28 21:51:48  vsnyder
! Make AT fields wider
!
! Revision 1.26  2017/11/28 01:46:49  vsnyder
! Don't compute Rel at Max Abs Diff = NaN for identical zero results
!
! Revision 1.25  2017/08/01 02:57:18  vsnyder
! Don't compare unequal-size records
!
! Revision 1.24  2017/08/01 02:32:36  vsnyder
! Use rmaxvg to decide whether to print summary
!
! Revision 1.23  2017/07/15 00:12:33  vsnyder
! Add -Z option
!
! Revision 1.22  2017/05/02 01:01:13  vsnyder
! Print maxdiff/maxval if -a
!
! Revision 1.21  2014/10/08 21:57:03  vsnyder
! Add 'RelMaxVal block' to summary header
!
! Revision 1.20  2014/10/08 20:22:21  vsnyder
! Add block label for RelMaxVal
!
! Revision 1.19  2013/08/06 23:14:31  vsnyder
! Remove dependence on machine module
!
! Revision 1.18  2009/04/13 20:43:17  pwagner
! Fixed a bug preventing macros file from using its own macros properly
!
! Revision 1.17  2007/06/26 00:33:43  vsnyder
! Print difference relative to max(maxval(abs(r1)),maxval(abs(r2))).
! Get rid of eps and outputs scaled by it.
!
! Revision 1.16  2007/06/19 00:29:39  vsnyder
! Print Max Abs Value if -a
!
! Revision 1.15  2007/06/08 22:46:03  vsnyder
! Revise output format
!
! Revision 1.14  2006/09/11 21:06:55  vsnyder
! In "all" format, print rel diff at max abs diff, abs diff at max rel diff
!
! Revision 1.13  2005/06/22 19:27:32  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.12  2004/09/23 23:33:39  vsnyder
! More futzing with -i option
!
! Revision 1.11  2004/09/23 23:01:43  vsnyder
! Don't print 'Identical' if file lengths different
!
! Revision 1.10  2004/09/17 20:59:34  vsnyder
! Add -i option
!
! Revision 1.9  2003/11/01 01:51:24  vsnyder
! Use single-precision epsilon for output computation
!
! Revision 1.8  2003/10/24 23:51:35  vsnyder
! Give AbsAtRmaxg and RelAtAmaxg initial NaN value
!
! Revision 1.7  2003/10/07 23:19:45  vsnyder
! Handle one-element dumps
!
! Revision 1.6  2003/10/07 02:04:03  vsnyder
! Show NaN relative differences differently
!
! Revision 1.5  2003/10/07 02:01:00  vsnyder
! Add option to print the version
!
! Revision 1.4  2003/10/07 01:57:17  vsnyder
! Detect NaN in relative difference
!
! Revision 1.3  2003/09/27 02:15:29  vsnyder
! Add computation of average and std dev
!
! Revision 1.2  2003/09/26 19:07:23  vsnyder
! Widen some formats
!
! Revision 1.1  2003/07/03 18:08:34  vsnyder
! Initial commit
!
