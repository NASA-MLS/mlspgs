program testRMS
   use DUMP_0, only: DUMP
   use MLSMessageModule, only: MLSMessageConfig
   use MLSStats1, only: STAT_T, &
     & ALLSTATS, DUMPSTAT=>DUMP, MLSMIN, MLSMAX, MLSMEAN, MLSSTDDEV, MLSRMS, STATISTICS
   use output_m, only: blanks, newline, output
  implicit none
  integer, parameter :: N=10
  real, dimension(N, N) :: array
  real :: min, max, mean, stddev, rms
  double precision, dimension(N, N) :: dblarray
  double precision :: dblmin, dblmax, dblmean, dblstddev, dblrms, eta
  type(STAT_T) :: statistic
  double precision, parameter :: X1 = -1.d0
  double precision, parameter :: X2 = 1.d0
  integer, parameter :: NCELLS = 12
  integer :: cell, status
  double precision, dimension(NCELLS-2) :: d1darray
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  array = 1.
  call dump ( array, &
    & 'all ones', stats=.true., rms=.true. )
  call allstats(array, max=max, min=min, mean=mean, stddev=stddev, rms=rms)
  call dumpstats(max, min, mean, stddev, rms)
  array(:,1:5) = 0.
  call dump ( array, &
    & 'half ones', stats=.true., rms=.true. )
  call allstats(array, max=max, min=min, mean=mean, stddev=stddev, rms=rms)
  call dumpstats(max, min, mean, stddev, rms)
  array(:,1:5) = -1.
  call dump ( array, &
    & 'half +/-', stats=.true., rms=.true. )
  call allstats(array, max=max, min=min, mean=mean, stddev=stddev, rms=rms)
  call dumpstats(max, min, mean, stddev, rms)
  array(:,6:10) = -999.99
  call dump ( array, &
    & 'half +/-', stats=.true., rms=.true., fillValue=-999.99 )
  call allstats(array, max=max, min=min, mean=mean, stddev=stddev, rms=rms, &
    & fillValue=-999.99)
  call dumpstats(max, min, mean, stddev, rms)
  dblarray = array
  dblarray(:,6:10) = -999.99d0
  call dump ( real(dblarray), &
    & 'half +/-', stats=.true., rms=.true., fillValue=-999.99 )
  call allstats(dblarray, max=dblmax, min=dblmin, mean=dblmean, &
    & stddev=dblstddev, rms=dblrms, fillValue=-999.99d0)
  call dumpstats(real(dblmax), real(dblmin), real(dblmean), real(dblstddev), &
    & real(dblrms))
  call output('mlsmin: ')
  call output(mlsmin(dblarray, fillValue=-999.99d0), advance='yes')
  call output('mlsmax: ')
  call output(mlsmax(dblarray, fillValue=-999.99d0), advance='yes')
  call output('mlsmean: ')
  call output(mlsmean(dblarray, fillValue=-999.99d0), advance='yes')
  call output('mlsstddev: ')
  call output(mlsstddev(dblarray, fillValue=-999.99d0), advance='yes')
  call output('mlsrms: ')
  call output(mlsrms(dblarray, fillValue=-999.99d0), advance='yes')
  ! Now let's do it all with user-defined type
  statistic%nbins = NCELLS
  allocate(statistic%bincount(NCELLS), stat=status)
  statistic%bounds = (/X1, X2/)
  do cell=1, NCELLS-2
    eta = (cell - 0.5)/ (NCELLS-2)
    d1darray(cell) = eta*X2 + (1-eta)*X1
    call statistics(d1darray(cell:cell), statistic)
  enddo
  call dumpstat(statistic)
  call statistics(d1darray, statistic)
  call dumpstat(statistic)
  do cell=1, NCELLS-2
    call statistics(d1darray, statistic)
  enddo
  call dumpstat(statistic)
contains
  subroutine dumpstats(max, min, mean, stddev, rms)
    ! dump the statistics
    ! Args
    real, intent(in)  :: max
    real, intent(in)  :: min
    real, intent(in)  :: mean
    real, intent(in)  :: stddev
    real, intent(in)  :: rms
    ! Executable
    call output('max: ')
    call output(max)
    call blanks(4)
    call output('min: ')
    call output(min)
    call blanks(4)
    call output('mean: ')
    call output(mean)
    call newline
    call output('stddev: ')
    call output(stddev)
    call blanks(4)
    call output('rms: ')
    call output(rms)
    call newline
  end subroutine dumpstats
end program testRMS
