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
    module procedure DUMP_3D_INTEGER
    module procedure DUMP_1D_CHAR, DUMP_2D_CHAR, DUMP_3D_CHAR
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  character, parameter :: AfterSub = '#'
  logical, parameter ::   DEEBUG = .false.

contains

  ! --------------------------------------------  DUMP_1D_CHAR  -----
  subroutine DUMP_1D_CHAR ( ARRAY, NAME, CLEAN )
    character(len=*), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    integer :: J, K
    logical :: MyClean
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do j = 1, size(array), 10
        if (.not. myClean) then
          if ( any(array(j:min(j+9, size(array))) /= ' ') ) then
            if ( numZeroRows /= 0 ) then
              call output ( j-1, places=max(4,ilog10(size(array))+1) )
              call output ( afterSub )
              call output ( ' ' )
              call output ( numZeroRows )
              call output ( ' rows of blanks not printed.', advance='yes' )
              numZeroRows = 0
            end if
            call output ( j, places=max(4,ilog10(size(array))+1) )
            call output ( afterSub )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+9, size(array))) /= ' ') ) then
          do k = j, min(j+9, size(array))
              call output ( array(k) // ' ' )
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      if ( numZeroRows /= 0 ) then
        call output ( j-1, places=max(4,ilog10(size(array))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of blanks not printed.', advance='yes' )
        numZeroRows = 0
      end if
    end if
  end subroutine DUMP_1D_CHAR

  ! --------------------------------------------  DUMP_2D_CHAR  -----
  subroutine DUMP_2D_CHAR ( ARRAY, NAME, CLEAN )
    character(len=*), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    integer :: I, J, K
    logical :: MyClean
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          if (.not. myClean) then
            if ( any(array(i,j:min(j+9, size(array,2))) /= ' ') ) then
              if ( numZeroRows /= 0 ) then
                call output ( i, places=max(4,ilog10(size(array,1))+1) )
                call output ( j-1, places=max(4,ilog10(size(array))+1) )
                call output ( afterSub )
                call output ( ' ' )
                call output ( numZeroRows )
                call output ( ' rows of blanks not printed.', advance='yes' )
                numZeroRows = 0
              end if
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+9, size(array,2))) /= ' ') ) then
            do k = j, min(j+9, size(array,2))
                call output ( array(i,k) // ' ' )
            end do
            call output ( '', advance='yes' )
          end if
        end do ! j
      end do ! i
      if ( numZeroRows /= 0 ) then
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
        call output ( j-1, places=max(4,ilog10(size(array))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of blanks not printed.', advance='yes' )
        numZeroRows = 0
      end if
    end if
  end subroutine DUMP_2D_CHAR

  ! ---------------------------------------------  DUMP_3D_CHAR  -----
  subroutine DUMP_3D_CHAR ( ARRAY, NAME, CLEAN )
    character(len=*), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many

    myClean = .false.
    if ( present(clean) ) myClean = clean
    call which_ints_are_it( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1,1), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, clean=clean )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, clean=clean )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 10
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+9, size(array,3))) /= ' ') ) then
                if ( numZeroRows /= 0 ) then
                  call output ( i, places=max(4,ilog10(size(array,1))+1) )
                  call output ( j, places=max(4,ilog10(size(array,2))+1) )
                  call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
                  call output ( afterSub )
                  call output ( ' ' )
                  call output ( numZeroRows )
                  call output ( ' rows of blanks not printed.', advance='yes' )
                  numZeroRows = 0
                end if
                call output ( i, max(4,ilog10(size(array,1))+1) )
                call output ( j, max(4,ilog10(size(array,2))+1) )
                call output ( k, max(4,ilog10(size(array,3))+1) )
                call output ( afterSub )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+9, size(array,3))) /= ' ') ) then
              do l = k, min(k+9, size(array,3))
                  call output ( array(i,j,l) // ' ' )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      if ( numZeroRows /= 0 ) then
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
        call output ( j-1, places=max(4,ilog10(size(array,2))+1) )
        call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of blanks not printed.', advance='yes' )
        numZeroRows = 0
      end if
    end if
  end subroutine DUMP_3D_CHAR

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
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1), '(1x,1pg13.6)', advance='yes' )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do j = 1, size(array), 5
        if (.not. myClean) then
          call output ( j, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+4, size(array))
          call output ( array(k), '(1x,1pg13.6)' )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    integer :: J, K
    logical :: MyClean
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do j = 1, size(array), 10
        if (.not. myClean) then
          if ( any(array(j:min(j+9, size(array))) /= 0) ) then
            if ( numZeroRows /= 0 ) then
              call output ( j-1, places=max(4,ilog10(size(array))+1) )
              call output ( afterSub )
              call output ( ' ' )
              call output ( numZeroRows )
              call output ( ' rows of zeros not printed.', advance='yes' )
              numZeroRows = 0
            end if
            call output ( j, places=max(4,ilog10(size(array))+1) )
            call output ( afterSub )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+9, size(array))) /= 0) ) then
          do k = j, min(j+9, size(array))
            if ( present(format) ) then
              call output ( array(k), format=format )
            else
              call output ( array(k), places=6 )
            end if
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      if ( numZeroRows /= 0 ) then
        call output ( j-1, places=max(4,ilog10(size(array))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of zeros not printed.', advance='yes' )
        numZeroRows = 0
      end if
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
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1), advance='yes' )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
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

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
  subroutine DUMP_2D_DOUBLE ( ARRAY, NAME, BLASE, CLEAN )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: BLASE
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K
    integer :: NumZeroRows
    double precision :: myBlase

    myBlase = 0.d0
    if ( present(Blase) ) myBlase = Blase

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else 
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), 5
            if (.not. myClean) then
              !call output ( i, max(4,ilog10(size(array,1))+1) )
              !call output ( j, max(4,ilog10(size(array,2))+1) )
              !call output ( afterSub )
              if ( any(array(i,j:min(j+4, size(array,2))) /= myBlase) ) then
                if ( numZeroRows /= 0 ) then
                  call output ( i, places=max(4,ilog10(size(array,1))+1) )
                  call output ( j-1, places=max(4,ilog10(size(array))+1) )
                  call output ( afterSub )
                  call output ( ' ' )
                  call output ( numZeroRows )
                  call output ( ' rows of ')
                  call output ( myBlase , advance='no' )
                  call output ( ' not printed.', advance='yes' )
                  numZeroRows = 0
                end if
                call output ( i, places=max(4,ilog10(size(array,1))+1) )
                call output ( j, places=max(4,ilog10(size(array,2))+1) )
                call output ( afterSub )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j:min(j+4, size(array,2))) /= myBlase) ) then
              do k = j, min(j+4, size(array,2))
                call output ( array(i,k), '(1x,1pg13.6)' )
              end do
              call output ( '', advance='yes' )
            endif
          end do
        end do
        if ( numZeroRows /= 0 ) then
          call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
          call output ( j-1, places=max(4,ilog10(size(array))+1) )
          call output ( afterSub )
          call output ( ' ' )
          call output ( numZeroRows )
          call output ( ' rows of ')                            
          call output ( myBlase , advance='no' )                
          call output ( ' not printed.', advance='yes' )        
          numZeroRows = 0
        end if
      else ! Dump the transpose
        call output ( ', transposed', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), 5
            !call output ( i, max(4,ilog10(size(array,1))+1) )
            !call output ( j, max(4,ilog10(size(array,2))+1) )
            !call output ( afterSub )
            if ( any(array(i:min(i+4, size(array,1)),j) /= myBlase) ) then  
              if ( numZeroRows /= 0 ) then                                  
                call output ( i, places=max(4,ilog10(size(array,1))+1) )    
                call output ( j-1, places=max(4,ilog10(size(array))+1) )    
                call output ( afterSub )                                    
                call output ( ' ' )                                         
                call output ( numZeroRows )                                 
                call output ( ' rows of ')                                  
                call output ( myBlase , advance='no' )                      
                call output ( ' not printed.', advance='yes' )              
                numZeroRows = 0                                             
              end if                                                        
              call output ( i, places=max(4,ilog10(size(array,1))+1) )      
              call output ( j, places=max(4,ilog10(size(array,2))+1) )      
              call output ( afterSub )                                      
            else                                                            
              numZeroRows = numZeroRows + 1                                 
            end if                                                          
            if ( myClean .or. any(array(i:min(i+4, size(array,1)),j) /= myBlase) ) then
              do k = i, min(i+4, size(array,1))
                call output ( array(k,j), '(1x,1pg13.6)' )
              end do
              call output ( '', advance='yes' )
            end if                                                          
          end do
        end do
      end if
      if ( numZeroRows /= 0 ) then                                  
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )  
        call output ( j-1, places=max(4,ilog10(size(array))+1) )    
        call output ( afterSub )                                    
        call output ( ' ' )                                         
        call output ( numZeroRows )                                 
        call output ( ' rows of ')                                  
        call output ( myBlase , advance='no' )                      
        call output ( ' not printed.', advance='yes' )              
        numZeroRows = 0                                             
      end if                                                        
    end if
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  subroutine DUMP_2D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    integer :: I, J, K
    logical :: MyClean
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          if (.not. myClean) then
            if ( any(array(i,j:min(j+9, size(array,2))) /= 0) ) then
              if ( numZeroRows /= 0 ) then
                call output ( i, places=max(4,ilog10(size(array,1))+1) )
                call output ( j-1, places=max(4,ilog10(size(array))+1) )
                call output ( afterSub )
                call output ( ' ' )
                call output ( numZeroRows )
                call output ( ' rows of zeros not printed.', advance='yes' )
                numZeroRows = 0
              end if
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+9, size(array,2))) /= 0) ) then
            do k = j, min(j+9, size(array,2))
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
      if ( numZeroRows /= 0 ) then
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
        call output ( j-1, places=max(4,ilog10(size(array))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of zeros not printed.', advance='yes' )
        numZeroRows = 0
      end if
    end if
  end subroutine DUMP_2D_INTEGER

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME, BLASE, CLEAN )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: BLASE
    logical, intent(in), optional :: CLEAN

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    double precision :: myBlase

    myBlase = 0.d0
    if ( present(Blase) ) myBlase = Blase
    myClean = .false.
    if ( present(clean) ) myClean = clean

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1,1), '(1x,1pg13.6)', advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name, clean=clean )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, clean=clean )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              !call output ( i, max(4,ilog10(size(array,1))+1) )
              !call output ( j, max(4,ilog10(size(array,2))+1) )
              !call output ( k, max(4,ilog10(size(array,3))+1) )
              !call output ( afterSub )
              if ( any(array(i,j,k:min(k+4, size(array,3))) /= myBlase) ) then
                if ( numZeroRows /= 0 ) then
                  call output ( i, places=max(4,ilog10(size(array,1))+1) )
                  call output ( j, places=max(4,ilog10(size(array,2))+1) )
                  call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
                  call output ( afterSub )
                  call output ( ' ' )
                  call output ( numZeroRows )
                  call output ( ' rows of ')
                  call output ( myBlase , advance='no' )
                  call output ( ' not printed.', advance='yes' )
                  numZeroRows = 0
                end if
                call output ( i, max(4,ilog10(size(array,1))+1) )
                call output ( j, max(4,ilog10(size(array,2))+1) )
                call output ( k, max(4,ilog10(size(array,3))+1) )
                call output ( afterSub )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+4, size(array,3))) /= myBlase) ) then
              do l = k, min(k+4, size(array,3))
                call output ( array(i,j,l), '(1x,1pg13.6)' )
              end do
              call output ( '', advance='yes' )
            endif
          end do
        end do
      end do
       if ( numZeroRows /= 0 ) then
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
        call output ( j-1, places=max(4,ilog10(size(array,2))+1) )
        call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of ')                              
        call output ( myBlase , advance='no' )                  
        call output ( ' not printed.', advance='yes' )          
        numZeroRows = 0
      end if
   end if
  end subroutine DUMP_3D_DOUBLE

  ! ---------------------------------------------  DUMP_3D_INTEGER  -----
  subroutine DUMP_3D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT )
    integer, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many

    myClean = .false.
    if ( present(clean) ) myClean = clean
    call which_ints_are_it( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      if ( present(name) ) then
        call output ( name )
        call output ( ' is ' )
      end if
      call output ( 'empty', advance='yes' )
    else if ( size(array) == 1 ) then
      if ( present(name) ) then
        call output ( name )
        if ( myClean ) call output ( ' \ 1 ' )
        call output ( ' ' )
      end if
      call output ( array(1,1,1), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, clean=clean )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, clean=clean )
    else
      if ( present(name) ) then 
        call output ( name )
        if ( myClean ) then 
          call output ( ' \ ' )
          call output ( size(array) )
        end if
        call output ( '', advance='yes' )
      end if
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 10
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+9, size(array,3))) /= 0) ) then
                if ( numZeroRows /= 0 ) then
                  call output ( i, places=max(4,ilog10(size(array,1))+1) )
                  call output ( j, places=max(4,ilog10(size(array,2))+1) )
                  call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
                  call output ( afterSub )
                  call output ( ' ' )
                  call output ( numZeroRows )
                  call output ( ' rows of zeros not printed.', advance='yes' )
                  numZeroRows = 0
                end if
                call output ( i, max(4,ilog10(size(array,1))+1) )
                call output ( j, max(4,ilog10(size(array,2))+1) )
                call output ( k, max(4,ilog10(size(array,3))+1) )
                call output ( afterSub )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+9, size(array,3))) /= 0) ) then
              do l = k, min(k+9, size(array,3))
                if ( present(format) ) then
                  call output ( array(i,j,l), format=format )
                else
                  call output ( array(i,j,l), places=6 )
                endif
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      if ( numZeroRows /= 0 ) then
        call output ( i-1, places=max(4,ilog10(size(array,1))+1) )
        call output ( j-1, places=max(4,ilog10(size(array,2))+1) )
        call output ( k-1, places=max(4,ilog10(size(array,3))+1) )
        call output ( afterSub )
        call output ( ' ' )
        call output ( numZeroRows )
        call output ( ' rows of zeros not printed.', advance='yes' )
        numZeroRows = 0
      end if
    end if
  end subroutine DUMP_3D_INTEGER

  integer function ilog10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ilog10

  ! Along with a (yet unwritten) fraternal subr. 'which_strings_are_it'
  ! this is possibly better aligned with MLSNumerics or MLSStrings
  ! themes; however, for now we leave it here
  ! -------------------------  which_ints_are_it  -----
  subroutine which_ints_are_it(ints, it, which, how_many, re_mainder, which_not)
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
      endif
      return
    endif
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
      endif
    enddo
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
        endif
        if ( present(which_not) ) then
          call output('which_not: ')
          call output(which_not, advance='yes')
        endif
        if ( present(re_mainder) ) then
          call output('re_mainder: ')
          call output(re_mainder, advance='yes')
        endif
    endif

  end subroutine which_ints_are_it

end module DUMP_0

! $Log$
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
