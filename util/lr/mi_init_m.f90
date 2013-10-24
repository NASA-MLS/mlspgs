module MI_INIT_M

! Machine-independent initialization

  use IO, only : INUNIT, PRUNIT, TBUNIT
  use MACHINE, only : IO_ERROR

  implicit NONE
  private

  public :: MI_INIT, OPEN_INPUT, OPEN_LISTING, OPEN_OUTPUT

contains

  subroutine MI_INIT
  ! Initialize various data structures
  end subroutine MI_INIT

  subroutine OPEN_INPUT ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )
    if ( opened ) return
    do
      if ( file /= ' ' ) then
        open ( inunit, file=file, status='OLD', form='FORMATTED', &
               iostat=IOSTAT )
        if ( iostat == 0 ) exit
        call io_error ( 'Unable to open input file', iostat, file )
      end if
      do
        write ( *, * ) 'Enter input file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do
    return
  end subroutine OPEN_INPUT

  subroutine OPEN_LISTING ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )

    do
      if ( file /= ' ' ) then
        open ( prunit, file=file, status='UNKNOWN', form='FORMATTED', &
               iostat=IOSTAT )
        if ( iostat == 0 ) exit
        call io_error ( 'Unable to open listing file', iostat, file )
      end if
      do
        write ( *, * ) 'Enter listing file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do
    return
  end subroutine OPEN_LISTING


  subroutine OPEN_OUTPUT ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )

    do
      if ( file /= ' ' ) then
        open ( tbunit, file=file, status='UNKNOWN', form='FORMATTED', &
               iostat=IOSTAT )
        if ( iostat == 0 ) exit
        call io_error ( 'Unable to open output file', iostat, file )
      end if
      do
        write ( *, * ) 'Enter output file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do
    return
  end subroutine OPEN_OUTPUT

end module MI_INIT_M
