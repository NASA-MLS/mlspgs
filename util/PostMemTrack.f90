! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program PostMemTrack

! Post-process the output created by the --memTrack option lf mlsl2.

! Try to match allocates and deallocates according to their names and sizes.

  type :: Record_t
    double precision :: Size = -1.0d0
    character(127) :: Where
  end type Record_t

  type(record_t), pointer :: Allocs(:), Temp(:)

  integer, parameter :: StartSize = 10000 ! of Allocs

  character(255) :: Line
  integer :: I, J, numAllocs = 0
  double precision :: Amount ! of an allocation or deallocation
  character(127) :: Where ! of an allocation or deallocation

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

  allocate ( allocs(startSize) )

o:do
    read ( *, '(a)', end=9 ) line
    if ( line(1:10) /= 'Tracking: ' ) cycle
    if ( line(11:20) == 'Allocated ' ) then
      do i = 1, numAllocs
        if ( allocs(i)%size < 0.0d0 ) exit
      end do
      if ( i > size(allocs) ) then
        temp => allocs
        allocate ( allocs(2*size(temp)) )
        allocs(:size(temp)) = temp
        deallocate ( temp )
      end if
      call Parse ( line(21:) )
      allocs(i)%size = amount
      allocs(i)%where = where
      numAllocs = max(numAllocs,i)
    else if ( line(11:22) == 'Deallocated ' ) then
      call parse ( line(23:) )
      do i = 1, numAllocs
        if ( allocs(i)%size == amount .and. &
          lowerCase(allocs(i)%where) == lowerCase(where) ) then
          allocs(i)%size = -1.0d0
          do j = numAllocs-1, 1, -1
            if ( allocs(j)%size >= 0.0d0 ) exit
          end do
          numAllocs = j
          cycle o
        end if
      end do
      write ( *, '(2a)' ) "***** Couldn't find an allocation for ", trim(line)
    else
      write ( *, '(2a)' ) '***** What is this: ', trim(line)
    end if
  end do o

9 continue
  do i = 1, size(allocs)
    if ( allocs(i)%size >= 0.0d0 ) &
      & write ( *, '("Leaked ", g16.8, " from ", a )' ) &
        & allocs(i)%size, trim(allocs(i)%where)
  end do

contains

  subroutine Parse ( Line )
  ! Get Size and Where from the Line
    character(len=*), intent(in) :: Line
    integer :: I1, I2
    character :: Scale ! First character of bytes, kb, mb, gb, tb
    ! Get Size
    read ( line, * ) amount, scale
    select case ( scale )
    case ( 'b', 'B' )
    case ( 'k', 'K' )
      amount = amount * 1024.0d0
    case ( 'm', 'M' )
      amount = amount * 1024.0d0**2
    case ( 'g', 'G' )
      amount = amount * 1024.0d0**3
    case ( 't', 'T' )
      amount = amount * 1024.0d0**4
    case default
      write ( *, '(2a)' ) '***** What scale is this: ', trim(line(8:))
    end select
    ! Get Where
    i1 = index(line,' for ')
    i2 = index(line,' total ', back=.true. )
    where = line(i1+5:i2-1)
  end subroutine Parse

  function LowerCase ( In )
    character(len=*), intent(in) :: In
    character(len=len(in)) :: LowerCase
    integer :: I
    lowerCase = in
    do i = 1, len(in)
      if ( in(i:i) >= 'A' .and. in(i:i) <= 'Z' ) &
        & lowerCase(i:i) = achar(iachar(in(i:i)) - (iachar('A') - iachar('a')))
    end do
  end function LowerCase

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end program PostMemTrack

! $Log$
