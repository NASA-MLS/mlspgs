! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program drsrang

  ! Demonstration driver for srang

  use MLSRandomNumber, only: MLS_RANDOM_SEED, SRANG
  use MLSStrings, only: LowerCASE

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
!---------------------------------------------------------------------------

  logical, parameter :: RANDOMRESTARTS = .false.
  logical, parameter :: READMYSEED = .true.
  logical, parameter :: ASKTOUSEINTRINSIC = .true.
  character(LEN=1)   :: answer
  integer :: seed(2), count, values(8)
  integer, parameter :: NCELLS = 12+2
  integer, parameter :: N = 10000
  real, parameter :: zero = 0.0
  real, parameter :: y1=-3.0, y2=3.0
  real :: stats(5), ytab(1)
  integer :: ihist(ncells)
!
if ( ASKTOUSEINTRINSIC ) then
   print*, 'Do you want mls_random_number switch to using f95 intrinsic?(y/n)'
   print*, '(Otherwise, I will use MATH77-derived methods)'
   read *, answer
   if ( LowerCase(answer) == 'y' ) then
      call MLS_RANDOM_SEED( MATH77_ranpack=.false. )
      print *, 'Switching to intrinsic random_number'
   endif
endif

if ( READMYSEED ) then
   print*, 'Please enter two integers to be used in forming seed'
   print*, '(if both 0, I will use default seed)'
   print*, '(if exactly one 0, I will let mls_random_seed choose seed)'
   print*, '(if neither 0, I will use your seed)'
   read *, seed
   if ( seed(1) /= 0 .and. seed(2) /= 0 ) then
      call MLS_RANDOM_SEED( pput=seed )
   elseif ( seed(1) == 0 .and. seed(2) == 0 ) then
      call MLS_RANDOM_SEED( gget=seed )
      print *, 'using default seed ', seed
!   elseif ( seed(1) /= 0 .or. seed(2) /= 0 ) then
   else
      call MLS_RANDOM_SEED( new_seed=seed )
      print *, 'new seed ', seed
   endif
endif
stats(1) = zero
do i=1, N
  if ( RANDOMRESTARTS .and. i > 1 .and. mod(i-1, 97)==0 ) &
    & call MLS_RANDOM_SEED()
  ytab(1) = srang()
  call sstat1(ytab(1), 1, stats, ihist, ncells, y1, y2)
enddo

call date_and_time(values=values)
call system_clock(count=count)
print *, 'date and time millisec ', values(8)
print *, 'system_clock count ', count
print *, 'Enforcing random restarts from same seed? ', RANDOMRESTARTS
print *, 'Gaussian random numbers from srang()'
call sstat2(stats, ihist, ncells, y1, y2)  
end program drsrang

! $Log$
