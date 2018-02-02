! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Test_Secant_Hunt

  use Secant_Hunt_m, only: Secant_Hunt

  integer :: Array(0:100)
  integer :: Failures = 0
  logical :: First = .false.         ! Find first of several equals
  integer :: Loc
  integer :: N
  integer :: Probes
  integer :: Target
  logical, parameter :: OK = .false. ! Show successes

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
!---------------------------------------------------------------------------

  array(0) = - huge(0)
  do
    read ( *, *, end=9 ) n, first
    read ( *, *, end=9 ), array(1:n)
    do target = array(1)-1, array(n)+1
      call secant_hunt ( array(1:n), target, loc, probes, first )
      if ( loc <= 0 ) then
        if ( target < array(1) ) then
          if ( OK ) print '(9(a,i0))', 'Target ', target, ' < Array(1) = ', array(1)
        else
          print '(9(a,i0))', 'FAIL LOC == 0 and Target ', target, ' >= Array(1)'
          failures = failures + 1
        end if
      else if ( loc == n .and. target >= array(loc) ) then
        if ( OK ) print '(9(a,i0))', 'Array(', loc, ') = ', array(loc), ' <= Target ', &
          & target, ' using ', probes, ' probes'
      else if ( loc == n .and. target < array(loc) ) then
        print '(9(a,i0))', 'FAIL Target ', target, ' < Array(', loc, ')'
        failures = failures + 1
      else if ( loc > n ) then
        print '(9(a,i0))', 'FAIL LOC ', loc, ' > N ', n
        failures = failures + 1
      else if ( array(loc) <= target .and. target <= array(loc+1) ) then
        if ( first .and. loc > 1 .and. target == array(loc-1) ) then
          print '(9(a,i0))', 'FAIL FIRST Array(', loc-1, ') = Array(', loc, ')'
          failures = failures + 1
        else if ( OK ) then
           print '(9(a,i0))', 'Array(', loc, ') = ', array(loc), ' <= Target ', target, &
          & ' < Array(', loc+1, ') = ', array(loc+1), ' using ', probes, ' probes'
        end if
      else if ( target < array(loc) ) then
        print '(9(a,i0))', 'FAIL low: Target ', target, ' < Array(', loc, ') = ', &
          & array(loc), ' using ', probes, ' probes'
        failures = failures + 1
      else if ( target > array(loc+1) ) then
        print '(9(a,i0))', 'FAIL high: Target ', target, ' > Array(', loc, ') = ', &
          & array(loc), ' using ', probes, ' probes'
        failures = failures + 1
      else
        print '(9(a,i0))', 'What happened?  LOC = ', loc, ' using ', probes, ' probes'
        failures = failures + 1
      end if
    end do
    if ( failures == 0 ) then
      print '(a)', 'No failures'
    else
      print '(i0,a)', failures, ' failures'
    end if

    failures = 0
    do target = array(1)-1, array(n)+1
      call secant_hunt ( real(array(1:n)), real(target), loc, probes, first )
      if ( loc <= 0 ) then
        if ( target < array(1) ) then
          if ( OK ) print '(9(a,i0))', 'Target ', target, ' < Array(1) = ', array(1)
        else
          print '(9(a,i0))', 'FAIL LOC == 0 and Target ', target, ' >= Array(1)'
        failures = failures + 1
        end if
      else if ( loc == n .and. target >= array(loc) ) then
        if ( OK ) print '(9(a,i0))', 'Array(', loc, ') = ', array(loc), ' <= Target ', &
          & target, ' using ', probes, ' probes'
      else if ( loc == n .and. target < array(loc) ) then
        print '(9(a,i0))', 'FAIL Target ', target, ' < Array(', loc, ')'
        failures = failures + 1
      else if ( loc > n ) then
        print '(9(a,i0))', 'FAIL LOC ', loc, ' > N ', n
        failures = failures + 1
      else if ( array(loc) <= target .and. target <= array(loc+1) ) then
        if ( first .and. loc > 1 .and. target == array(loc-1) ) then
          print '(9(a,i0))', 'FAIL FIRST Array(', loc-1, ') = Array(', loc, ')'
          failures = failures + 1
        else if ( OK ) then
          print '(9(a,i0))', 'Array(', loc, ') = ', array(loc), ' <= Target ', target, &
            & ' < Array(', loc+1, ') = ', array(loc+1), ' using ', probes, ' probes'
        end if
      else if ( target < array(loc) ) then
        print '(9(a,i0))', 'FAIL low: Target ', target, ' < Array(', loc, ') = ', &
          & array(loc), ' using ', probes, ' probes'
        failures = failures + 1
      else if ( target > array(loc+1) ) then
        print '(9(a,i0))', 'FAIL high: Target ', target, ' > Array(', loc, ') = ', &
          & array(loc), ' using ', probes, ' probes'
        failures = failures + 1
      else
        print '(9(a,i0))', 'What happened?  LOC = ', loc, ' using ', probes, ' probes'
        failures = failures + 1
      end if
    end do
    if ( failures == 0 ) then
      print '(a)', 'No failures'
    else
      print '(i0,a)', failures, ' failures'
    end if
  end do

9 continue

end program Test_Secant_Hunt

! $Log$
