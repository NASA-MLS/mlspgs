! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DUMP_0

! Low-level dump routines -- for some arrays of intrinsic type.

  use OUTPUT_M, only: OUTPUT

  implicit NONE
  private
  public :: AfterSub, DUMP

  interface DUMP
    module procedure DUMP_1D_DOUBLE, DUMP_1D_INTEGER, DUMP_1D_LOGICAL
    module procedure DUMP_2D_DOUBLE, DUMP_2D_INTEGER, DUMP_3D_DOUBLE
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  character, parameter :: AfterSub = '#'

contains

  ! ---------------------------------------------  DUMP_1D_DOUBLE  -----
  subroutine DUMP_1D_DOUBLE ( ARRAY, NAME, CLEAN )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1), '(1x,1pg13.6)', advance='yes' )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do j = 1, size(array), 5
        if (.not. myClean) then
          call output ( j, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        endif
        do k = j, min(j+4, size(array))
          call output ( array(k), '(1x,1pg13.6)' )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME, CLEAN )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do j = 1, size(array), 10
        if (.not. myClean) then
          call output ( j, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        endif
        do k = j, min(j+9, size(array))
          call output ( array(k), 6 )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_INTEGER

  ! --------------------------------------------  DUMP_1D_LOGICAL ----
  subroutine DUMP_1D_LOGICAL ( ARRAY, NAME, CLEAN )
    logical, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do j = 1, size(array), 34
        if (.not. myClean) then
          call output ( j, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        endif
        do k = j, min(j+33, size(array))
          call output ( array(k) )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_LOGICAL

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
  subroutine DUMP_2D_DOUBLE ( ARRAY, NAME, CLEAN )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 5
          if (.not. myClean) then
            call output ( i, max(4,ilog10(size(array,1))+1) )
            call output ( j, max(4,ilog10(size(array,2))+1) )
            call output ( afterSub )
          endif
          do k = j, min(j+4, size(array,2))
            call output ( array(i,k), '(1x,1pg13.6)' )
          end do
          call output ( '', advance='yes' )
        end do
      end do
    end if
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  subroutine DUMP_2D_INTEGER ( ARRAY, NAME, CLEAN )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          if (.not. myClean) then
            call output ( i, max(4,ilog10(size(array,1))+1) )
            call output ( j, max(4,ilog10(size(array,2))+1) )
            call output ( afterSub )
          end if
          do k = j, min(j+9, size(array,2))
            call output ( array(i,k), 6 )
          end do
          call output ( '', advance='yes' )
        end do
      end do
    end if
  end subroutine DUMP_2D_INTEGER

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME, CLEAN )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K, L

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' ' )
      end if
      call output ( array(1,1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name )
    else
      if ( present(name) ) call output ( name, advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              call output ( i, max(4,ilog10(size(array,1))+1) )
              call output ( j, max(4,ilog10(size(array,2))+1) )
              call output ( k, max(4,ilog10(size(array,3))+1) )
              call output ( afterSub )
            endif
            do l = k, min(k+4, size(array,3))
              call output ( array(i,j,l), '(1x,1pg13.6)' )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end do
    end if
  end subroutine DUMP_3D_DOUBLE

  integer function ilog10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ilog10

end module DUMP_0

! $Log$
! Revision 2.6  2001/05/08 17:21:02  livesey
! Added a `clean' option to the array dumps.  This omits the indices at
! the start, making it easier for other programs to read output.
!
! Revision 2.5  2001/05/03 02:12:34  vsnyder
! Insert copyright notice, clean up CVS stuff, cosmetics
!
! Revision 2.4  2001/03/10 03:39:58  vsnyder
! Improve handling of "name" if size==1 or size==0
!
! Revision 2.3  2001/03/02 01:32:08  livesey
! Handles larger arrays better
!
! Revision 2.2  2001/02/28 21:35:27  livesey
! Added dump logical 1d
!
! Revision 2.1  2000/09/13 20:38:50  vsnyder
! Initial code
!
