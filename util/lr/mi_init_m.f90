! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MI_INIT_M

! Machine-independent initialization

  use IO, only : INUNIT, PRUNIT, TBUNIT

  implicit NONE
  private

  public :: MI_INIT, OPEN_INPUT, OPEN_LISTING, OPEN_OUTPUT

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine MI_INIT
  ! Initialize various data structures
  end subroutine MI_INIT

  subroutine OPEN_INPUT ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    character(len=127) :: MSG ! in case of I/O error
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )
    if ( opened ) return
    do
      if ( file /= ' ' ) then
        open ( inunit, file=file, status='OLD', form='FORMATTED', &
               iostat=IOSTAT, iomsg=MSG )
        if ( iostat == 0 ) exit
        write ( *, '(a,i0/3a/a)' ) 'Unable to open input file, IOSTAT = ', &
          & iostat, 'File: "', trim(file), '"', trim(msg)
      end if
      do
        write ( *, * ) 'Enter input file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do

  end subroutine OPEN_INPUT

  subroutine OPEN_LISTING ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    character(len=127) :: MSG ! in case of I/O error
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )

    do
      if ( file /= ' ' ) then
        open ( prunit, file=file, status='UNKNOWN', form='FORMATTED', &
               iostat=IOSTAT, iomsg=MSG )
        if ( iostat == 0 ) exit
        write ( *, '(a,i0/3a/a)' ) 'Unable to open listing file, IOSTAT = ', &
          & iostat, 'File: "', trim(file), '"', trim(msg)
      end if
      do
        write ( *, * ) 'Enter listing file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do

  end subroutine OPEN_LISTING


  subroutine OPEN_OUTPUT ( FILE )
    character(len=*), intent(inout) :: FILE

    integer :: IOSTAT
    character(len=127) :: MSG ! in case of I/O error
    logical :: OPENED

    inquire ( unit=inunit, opened=opened )

    do
      if ( file /= ' ' ) then
        open ( tbunit, file=file, status='UNKNOWN', form='FORMATTED', &
               iostat=IOSTAT, iomsg=MSG )
        if ( iostat == 0 ) exit
        write ( *, '(a,i0/3a/a)' ) 'Unable to open output file, IOSTAT = ', &
          & iostat, 'File: "', trim(file), '"', trim(msg)
      end if
      do
        write ( *, * ) 'Enter output file name: '
        read ( *, '(a)', iostat = iostat ) file
        if ( iostat /= 0 ) stop
        if ( file /= ' ' ) exit
        write ( *, * ) 'Nothing entered.'
      end do
    end do

  end subroutine OPEN_OUTPUT

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MI_INIT_M

! $Log$
