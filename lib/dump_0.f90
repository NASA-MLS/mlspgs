! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DUMP_0

! Low-level dump routines -- for some arrays of intrinsic type.

  use MLSStrings, only: GetStringElement
  use OUTPUT_M, only: BLANKS, OUTPUT

  implicit NONE
  private
  public :: AfterSub, DUMP, DUMP_NAME_V_PAIRS

  interface DUMP        ! dump n-d arrays of homogeneous type
    module procedure DUMP_1D_CHAR, DUMP_1D_COMPLEX, DUMP_1D_DCOMPLEX
    module procedure DUMP_1D_DOUBLE, DUMP_1D_INTEGER
    module procedure DUMP_1D_LOGICAL, DUMP_1D_REAL
    module procedure DUMP_2D_CHAR, DUMP_2D_COMPLEX, DUMP_2D_DCOMPLEX
    module procedure DUMP_2D_DOUBLE, DUMP_2D_INTEGER
    module procedure DUMP_2D_LOGICAL, DUMP_2D_REAL
    module procedure DUMP_3D_CHAR, DUMP_3D_DOUBLE, DUMP_3D_INTEGER
    module procedure DUMP_3D_REAL
  end interface
  interface DUMP_NAME_V_PAIRS   ! dump name-value pairs, names in string list
    module procedure DUMP_NAME_V_PAIRS_DOUBLE, DUMP_NAME_V_PAIRS_INTEGER
    module procedure DUMP_NAME_V_PAIRS_REAL
  end interface

  interface SAY_FILL
    module procedure SAY_FILL_CHAR, SAY_FILL_DOUBLE, SAY_FILL_INT, SAY_FILL_REAL
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  character, parameter :: AfterSub = '#'
  logical, parameter ::   DEEBUG = .false.

  character(*), parameter :: MyFormatDefault = '(1pg14.6)'
  character(*), parameter :: MyFormatDefaultCmplx = &
    & '(1x,"(",1pg13.6,",",1pg13.6,")")'

contains

  ! -----------------------------------------------  DUMP_1D_CHAR  -----
  subroutine DUMP_1D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN )
    character(len=*), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN

    integer :: J, K
    logical :: MyClean
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), 10
        if (.not. myClean) then
          if ( any(array(j:min(j+9, size(array))) /= myFillValue) ) then
            call say_fill ( (/ j-1, size(array) /), numZeroRows, myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+9, size(array))) /= myFillValue) ) then
          do k = j, min(j+9, size(array))
              call output ( array(k) // ' ' )
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      call say_fill ( (/ j-10, size(array) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_1D_CHAR

  ! --------------------------------------------  DUMP_1D_COMPLEX  -----
  subroutine DUMP_1D_COMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: MyClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 3
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_COMPLEX

  ! -------------------------------------------  DUMP_1D_DCOMPLEX  -----
  subroutine DUMP_1D_DCOMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: MyClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 3
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DCOMPLEX

 ! ---------------------------------------------  DUMP_1D_DOUBLE  -----
  subroutine DUMP_1D_DOUBLE ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: MyClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 5
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT, WIDTH )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?

    integer :: J, K
    logical :: MyClean
    integer :: MyWidth
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          if ( any(array(j:min(j+myWidth-1, size(array))) /= 0) ) then
            call say_fill ( (/ j-1, size(array) /), numZeroRows, 0, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+myWidth-1, size(array))) /= 0) ) then
          do k = j, min(j+myWidth-1, size(array))
            if ( present(format) ) then
              call output ( array(k), format=format )
            else
              call output ( array(k), places=6 )
            end if
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      call say_fill ( (/ j-myWidth, size(array) /), numZeroRows, 0 )
    end if
  end subroutine DUMP_1D_INTEGER

  ! ----------------------------------------------  DUMP_1D_LOGICAL ----
  subroutine DUMP_1D_LOGICAL ( ARRAY, NAME, CLEAN )
    logical, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: J, K

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), 34
        if (.not. myClean) then
          call output ( j, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+33, size(array))
          call output ( array(k) )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_LOGICAL

  ! -----------------------------------------------  DUMP_1D_REAL  -----
  subroutine DUMP_1D_REAL ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 5
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_REAL

  ! -----------------------------------------------  DUMP_2D_CHAR  -----
  subroutine DUMP_2D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN )
    character(len=*), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN

    integer :: I, J, K
    logical :: MyClean
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          if (.not. myClean) then
            if ( any(array(i,j:min(j+9, size(array,2))) /= myFillValue) ) then
              call say_fill ( (/ i-1, size(array,1), j, size(array,2) /), &
                & numZeroRows, myFillValue, inc=1 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+9, size(array,2))) /= myFillValue) ) then
            do k = j, min(j+9, size(array,2))
                call output ( array(i,k) // ' ' )
            end do
            call output ( '', advance='yes' )
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size(array,1), j-10, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_CHAR

  ! --------------------------------------------  DUMP_2D_COMPLEX  -----
  subroutine DUMP_2D_COMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: myWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myWidth = 3
    if ( present(width) ) myWidth = width

    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), format=myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=myFormat )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), myWidth
            if (.not. myClean) then
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            end if
            do k = j, min(j+myWidth-1, size(array,2))
              call output ( array(i,k), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), width
            call output ( i, places=max(4,ilog10(size(array,1))+1) )      
            call output ( j, places=max(4,ilog10(size(array,2))+1) )      
            call output ( afterSub )                                      
            do k = i, min(i+width-1, size(array,1))
              call output ( array(k,j), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end if
    end if
  end subroutine DUMP_2D_COMPLEX

  ! --------------------------------------------  DUMP_2D_COMPLEX  -----
  subroutine DUMP_2D_DCOMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: myWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myWidth = 3
    if ( present(width) ) myWidth = width

    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), format=myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=myFormat )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), myWidth
            if (.not. myClean) then
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            end if
            do k = j, min(j+myWidth-1, size(array,2))
              call output ( array(i,k), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), width
            call output ( i, places=max(4,ilog10(size(array,1))+1) )      
            call output ( j, places=max(4,ilog10(size(array,2))+1) )      
            call output ( afterSub )                                      
            do k = i, min(i+width-1, size(array,1))
              call output ( array(k,j), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end if
    end if
  end subroutine DUMP_2D_DCOMPLEX

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
  subroutine DUMP_2D_DOUBLE ( ARRAY, NAME, FILLVALUE, CLEAN, FORMAT )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myFillValue = 0.0d0
    if ( present(FillValue) ) myFillValue = FillValue

    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), 5
            if (.not. myClean) then
              if ( any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                  & numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
              do k = j, min(j+4, size(array,2))
                call output ( array(i,k), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
        call say_fill ( (/ i-1, size(array,1), j-5, size(array,2) /), &
          & numZeroRows, myFillValue )
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), 5
            if ( any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then  
              call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                & numZeroRows, myFillValue, inc=3 )
            else                                                            
              numZeroRows = numZeroRows + 1                                 
            end if                                                          
            if ( myClean .or. any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then
              do k = i, min(i+4, size(array,1))
                call output ( array(k,j), myFormat )
              end do
              call output ( '', advance='yes' )
            end if                                                          
          end do
        end do
      end if
      call say_fill ( (/ i-5, size(array,1), j-1, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  subroutine DUMP_2D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT, WIDTH )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?

    integer :: I, J, K
    logical :: MyClean
    integer :: MyWidth
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=format, width=width )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), myWidth
          if (.not. myClean) then
            if ( any(array(i,j:min(j+myWidth-1, size(array,2))) /= 0) ) then
              call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                & numZeroRows, 0, inc=3 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+myWidth-1, size(array,2))) /= 0) ) then
            do k = j, min(j+myWidth-1, size(array,2))
              if ( present(format) ) then
                call output ( array(i,k), format=format )
              else
                call output ( array(i,k), places=6 )
              end if
            end do
            call output ( '', advance='yes' )
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size(array,1), j-myWidth, size(array,2) /), &
        & numZeroRows, 0 )
    end if
  end subroutine DUMP_2D_INTEGER

  ! --------------------------------------------  DUMP_2D_LOGICAL  -----
  subroutine DUMP_2D_LOGICAL ( ARRAY, NAME, CLEAN )
    logical, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    integer :: I, J, K
    logical :: MyClean
    integer, parameter :: MyWidth = 34

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), myWidth
          if (.not. myClean) then
            call output ( i, places=max(4,ilog10(size(array,1))+1) )
            call output ( j, places=max(4,ilog10(size(array,2))+1) )
            call output ( afterSub )
          end if
          do k = j, min(j+myWidth-1, size(array,2))
            call output ( array(i,k) )
          end do
          call output ( '', advance='yes' )
        end do ! j
      end do ! i
    end if
  end subroutine DUMP_2D_LOGICAL

  ! -----------------------------------------------  DUMP_2D_REAL  -----
  subroutine DUMP_2D_REAL ( ARRAY, NAME, FILLVALUE, CLEAN, FORMAT )
    real, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: NumZeroRows
    real :: myFillValue
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myFillValue = 0.0e0
    if ( present(FillValue) ) myFillValue = FillValue

    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), 5
            if (.not. myClean) then
              if ( any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                  & numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
              do k = j, min(j+4, size(array,2))
                call output ( array(i,k), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
        call say_fill ( (/ i-1, size(array,1), j-5, size(array,2) /), &
          & numZeroRows, myFillValue )
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), 5
            if ( any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then  
              call say_fill ( (/ i-1, size(array,1), j, size(array,2) /), &
                & numZeroRows, myFillValue, inc=1 )
            else                                                            
              numZeroRows = numZeroRows + 1                                 
            end if                                                          
            if ( myClean .or. any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then
              do k = i, min(i+4, size(array,1))
                call output ( array(k,j), myFormat )
              end do
              call output ( '', advance='yes' )
            end if                                                          
          end do
        end do
      end if
      call say_fill ( (/ i-5, size(array,1), j-1, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_REAL

  ! -----------------------------------------------  DUMP_3D_CHAR  -----
  subroutine DUMP_3D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN )
    character(len=*), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean
    call which_ints_are_it( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, fillValue=fillValue, &
        & clean=clean )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 10
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+9, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+9, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+9, size(array,3))
                  call output ( array(i,j,l) // ' ' )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-10, size(array,3) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_3D_CHAR

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME, FILLVALUE, CLEAN, FORMAT )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: myFormat

    myFillValue = 0.d0
    if ( present(FillValue) ) myFillValue = FillValue
    myClean = .false.
    if ( present(clean) ) myClean = clean
    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name, clean=clean )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+4, size(array,3))
                call output ( array(i,j,l), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-5, size(array,3) /), numZeroRows, myFillValue )
   end if
  end subroutine DUMP_3D_DOUBLE

  ! --------------------------------------------  DUMP_3D_INTEGER  -----
  subroutine DUMP_3D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT, WIDTH )
    integer, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?

    integer :: I, J, K, L
    logical :: myClean
    integer :: MyWidth
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width
    call which_ints_are_it( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, clean=clean, &
      & format=format )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, clean=clean, format=format )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), myWidth
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+myWidth-1, size(array,3))) /= 0) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, 0, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+myWidth-1, size(array,3))) /= 0) ) then
              do l = k, min(k+myWidth-1, size(array,3))
                if ( present(format) ) then
                  call output ( array(i,j,l), format=format )
                else
                  call output ( array(i,j,l), places=6 )
                end if
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-myWidth, size(array,3) /), numZeroRows, 0 )
    end if
  end subroutine DUMP_3D_INTEGER

  ! ---------------------------------------------  DUMP_3D_REAL  -----
  subroutine DUMP_3D_REAL ( ARRAY, NAME, FILLVALUE, CLEAN, FORMAT )
    real, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    real    :: myFillValue
    character(len=64) :: MyFormat

    myFillValue = 0.e0
    if ( present(FillValue) ) myFillValue = FillValue
    myClean = .false.
    if ( present(clean) ) myClean = clean
    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name, clean=clean )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+4, size(array,3))
                call output ( array(i,j,l), myFormat )
              end do
              call output ( '', advance='yes' )
            endif
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-5, size(array,3) /), numZeroRows, myFillValue )
   end if
  end subroutine DUMP_3D_REAL


  ! -----------------------------------  DUMP_NAME_V_PAIRS_DOUBLE  -----
  subroutine DUMP_NAME_V_PAIRS_DOUBLE ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    double precision, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'double # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_DOUBLE

  ! ----------------------------------  DUMP_NAME_V_PAIRS_INTEGER  -----
  subroutine DUMP_NAME_V_PAIRS_INTEGER ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    integer, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'integer # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format=format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_INTEGER

  ! -------------------------------------  DUMP_NAME_V_PAIRS_REAL  -----
  subroutine DUMP_NAME_V_PAIRS_REAL ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    real, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'single # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_REAL

  ! ------------------------------------------------------  Empty  -----
  subroutine Empty ( Name )
    character(len=*), intent(in), optional :: Name

    if ( present(name) ) then
      call output ( name )
      call output ( ' is ' )
    end if
    call output ( 'empty', advance='yes' )

  end subroutine Empty

  ! -----------------------------------------------------  ILOG10  -----
  integer function ILOG10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ILOG10

  ! ----------------------------------------------  Name_And_Size  -----
  subroutine Name_And_Size ( Name, Clean, Size )
    character(len=*), intent(in), optional :: Name
    logical, intent(in) :: Clean
    integer, intent(in) :: Size

    if ( present(name) ) then
      call output ( name )
      if ( clean ) then 
        call output ( trim(" \ ") ) ! This goofiness is to outwit an incorrect
                                    ! Intel compiler.
        call output ( size )
      end if
      if ( size == 1 ) call output ( ' ' )
    end if

  end subroutine Name_And_Size

  ! ----------------------------------------------  Say_Fill_Char  -----
  subroutine Say_Fill_Char ( Subs, NumZeroRows, Fill, Inc  )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    character(len=*), intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows  )
      call output ( '"', advance='no' )
      call output ( trim(fill), advance='no' )
      call output ( '" not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Char

  ! --------------------------------------------  Say_Fill_Double  -----
  subroutine Say_Fill_Double ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    double precision, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Double

  ! -----------------------------------------------  Say_Fill_Int  -----
  subroutine Say_Fill_Int ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    integer, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Int

  ! ----------------------------------------------  Say_Fill_Real  -----
  subroutine Say_Fill_Real ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    real, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Real

  ! ----------------------------------------------------  Say_Subs -----
  subroutine Say_Subs ( Subs, NumZeroRows )
    integer, intent(in) :: Subs(:)
    integer, intent(in) :: NumZeroRows
    call say_subs_only ( subs )
    call output ( ' ' )
    call output ( numZeroRows )
    call output ( ' rows of ', advance='no' )
  end subroutine Say_Subs

  ! -----------------------------------------------  Say_Subs_Only -----
  subroutine Say_Subs_Only ( Subs )
    integer, intent(in) :: Subs(:)
    integer :: I
    do i = 1, size(subs), 2
      call output ( subs(i), places=max(4,ilog10(subs(i+1))+1) )
    end do
    call output ( afterSub )
  end subroutine Say_Subs_Only

  ! ------------------------------------------  WHICH_INTS_ARE_IT  -----
  ! Along with a (yet unwritten) fraternal subr. 'which_strings_are_it'
  ! this is possibly better aligned with MLSNumerics or MLSStrings
  ! themes; however, for now we leave it here
  subroutine WHICH_INTS_ARE_IT ( INTS, IT, WHICH, HOW_MANY, RE_MAINDER, WHICH_NOT )
    ! Return which i of ints[i] = it
    ! optionally, return also how many of them do
    ! which_not of them (which don't)
    ! and the re_mainder of the ints != it
    ! e.g. given ints = /(4, 3, 1, 2, 1, 3 )/ and it = 1
    ! produces which = /(3, 5)/, 
    !      which_not = /(1, 2, 4, 6)/, 
    !       how_many = 2,
    !     re_mainder = /(4, 3, 2, 3)/
    
    ! This may be useful in reshaping an array to suppress any dims 
    ! that are identically 1

    ! Formal arguments
    integer, intent(in), dimension(:)  ::           ints
    integer, intent(in)                ::           it
    integer, intent(out), dimension(:) ::           which
    integer, intent(out), optional ::               how_many
    integer, intent(out), dimension(:), optional :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    
    if ( size(ints) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      if ( present(re_mainder) ) re_mainder = 0
      if ( DEEBUG ) then
        call output('size of ints or which too small', advance='yes')
        call output('ints: ')
        call output(ints, advance='yes')
        call output('which: ')
        call output(which, advance='yes')
      end if
      return
    end if
    i_which = 0
    i_re_mainder = 0
    do i=1, size(ints)
      if ( ints(i) == it ) then
        i_which = i_which+1
        which(min(size(which), i_which)) = i
      else
        i_re_mainder = i_re_mainder+1
        if ( present(which_not) ) &
          & re_mainder(min(size(which_not), i_re_mainder)) = i
        if ( present(re_mainder) ) &
          & re_mainder(min(size(re_mainder), i_re_mainder)) = ints(i)
      end if
    end do
    if ( present(how_many) ) how_many = i_which
    if ( DEEBUG ) then
        call output('ints: ')
        call output(ints, advance='yes')
        call output('it: ')
        call output(it, advance='yes')
        call output('which: ')
        call output(which, advance='yes')
        if ( present(how_many) ) then
          call output('how_many: ')
          call output(how_many, advance='yes')
        end if
        if ( present(which_not) ) then
          call output('which_not: ')
          call output(which_not, advance='yes')
        end if
        if ( present(re_mainder) ) then
          call output('re_mainder: ')
          call output(re_mainder, advance='yes')
        end if
    end if

  end subroutine WHICH_INTS_ARE_IT

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DUMP_0

! $Log$
! Revision 2.30  2003/09/19 02:00:14  vsnyder
! More about the goofy Intel compiler
!
! Revision 2.29  2003/09/15 17:43:41  livesey
! Cosmetic change for fussy (and wrong) intel compiler
!
! Revision 2.28  2003/09/06 00:48:40  vsnyder
! Specify default formats with a module parameter instead of literals.
! Change default (1x,1pg13.6) to (1pg14.6) to avoid problems with length
! calculation in output_m.
!
! Revision 2.27  2003/08/08 20:45:42  vsnyder
! Made say_fill_* generic, made them test for numZeroRows, and made them
! optionally do say_subs_only.  This simplified several dump routines.
! Added optional FORMAT arguments in several more routines.
!
! Revision 2.26  2003/07/04 02:41:33  vsnyder
! Substantial simplification by putting little things into subroutines
!
! Revision 2.25  2003/07/02 01:07:27  vsnyder
! Add complex output
!
! Revision 2.24  2003/05/21 19:20:40  vsnyder
! Start a new line after \ 1 if size==1 and clean
!
! Revision 2.23  2003/05/06 00:15:03  pwagner
! Fixed incompatibility with FilterShapes
!
! Revision 2.20.2.4  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.20.2.3  2003/04/18 20:26:05  vsnyder
! Add WIDTH and FORMAT arguments to 1D_REAL and 1D_DOUBLE
!
! Revision 2.20.2.2  2003/03/27 23:18:33  vsnyder
! Put new-lines in better places
!
! Revision 2.20.2.1  2003/03/14 00:25:47  vsnyder
! Add Dump_2D_Logical, cosmetic changes
!
! Revision 2.20  2002/12/02 23:34:14  pwagner
! Now can dump name/value pairs
!
! Revision 2.19  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.17  2002/02/14 23:21:18  vsnyder
! Work on dumping masks
!
! Revision 2.16  2001/12/08 00:47:51  pwagner
! Added dump_1d_real for s.p. arrays
!
! Revision 2.15  2001/11/29 23:50:53  pwagner
! Added optional blase arg to dump_nd_char; fixed bug where optional
! format not passed from dump_3d_int
!
! Revision 2.14  2001/11/28 23:32:01  livesey
! Fixed bug where dump_2d_integer didn't pass format to 1d dump.
!
! Revision 2.13  2001/10/25 23:30:39  pwagner
! Improved dump_nd_double to skip rows (e.g., of zeros)
!
! Revision 2.12  2001/10/24 18:11:14  pwagner
! which_ints_are_it now works properly
!
! Revision 2.11  2001/10/23 22:40:37  pwagner
! Now dumps 1d,2d,3d char arrays and 3d ints
!
! Revision 2.10  2001/09/28 22:43:20  vsnyder
! Don't print rows of zeroes
!
! Revision 2.9  2001/09/11 22:52:32  livesey
! Added printing of sizes
!
! Revision 2.8  2001/05/11 22:44:54  vsnyder
! Print transpose of 2d-double if it would take fewer lines.  Get rid of
! double printing of "without mask"
!
! Revision 2.7  2001/05/08 20:27:24  vsnyder
! Added an optional 'format' argument in a few more places
!
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
