! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Sparse_m

! Representation and some operations on sparse matrices

  use MLSKinds, only: RP

  implicit NONE

  private
  public :: Sparse_t
  public :: Create_Sparse, Destroy_Sparse, Add_Element, Resize_Sparse
  public :: Row_Dot_Vec, Sparse_Dot_Vec, Vec_Dot_Col, Vec_Dot_Sparse
  public :: Dump_Sparse, Dump

  type :: Sparse_Element_t ! One element in a sparse matrix
    real(rp) :: V
    integer :: R   ! In which row is the element
    integer :: C   ! In which col is the element
    integer :: NR  ! Next element in same row
    integer :: NC  ! Next element in same col
  end type Sparse_Element_t

  type :: Sparse_t ! Representation for a sparse matrix
    integer :: NE = 0   ! Number of elements actually used, <= size(E)
    integer :: What = 0 ! A string index for dumps
    ! Rows and columns are circular lists, so last element points to first:
    integer, allocatable :: Rows(:) ! Last element in each row
    integer, allocatable :: Cols(:) ! Last element in each col
    integer, allocatable :: LBnd(:) ! Column sub-dimension lower bounds
    integer, allocatable :: UBnd(:) ! Column sub-dimension upper bounds
    type(sparse_element_t), allocatable :: E(:) ! nonzero elements
  contains
    procedure :: Create => Create_Sparse
    procedure :: Add_Element_Stru
    procedure :: Add_Element_Value
    procedure :: Add_Elements_Stru
    generic :: Add_Element => Add_Element_Stru, Add_Element_Value, Add_Elements_Stru
    procedure :: Resize => Resize_Sparse
    procedure :: Destroy => Destroy_Sparse
    procedure :: Dump => Dump_Sparse
    procedure :: Row_Dot_Vec
    procedure :: Sparse_Dot_Vec
    procedure, pass(sparse) :: Vec_Dot_Col
    procedure, pass(sparse) :: Vec_Dot_Sparse
  end type Sparse_t

  interface Add_Element
    module procedure Add_Element_Stru, Add_Element_Value, Add_Elements_Stru
  end interface Add_Element

  interface Dump
    module procedure Dump_Sparse
  end interface Dump

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Create_Sparse ( Sparse, NR, NC, NE, UBnd, LBnd, What )
    ! Create a sparse matrix representation
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    class(sparse_t), intent(inout) :: Sparse  ! The matrix
    integer, intent(in) :: NR                 ! Number of rows
    integer, intent(in) :: NC                 ! Number of cols
    integer, intent(in), optional :: NE       ! Initial size ( sparse%e )
    integer, intent(in), optional :: UBnd(:)  ! Column sub-dimensions upper bounds
    integer, intent(in), optional :: LBnd(:)  ! Column sub-dimensions lower bounds
                                              ! Not used if UBnd not present
    integer, intent(in), optional :: What     ! A string index, for dumps

    character(127) :: Msg
    integer :: Stat

    call destroy_sparse ( sparse )
    call allocate_test ( sparse%rows, nr, "Sparse%Rows", moduleName, fill=0 )
    call allocate_test ( sparse%cols, nc, "Sparse%Cols", moduleName, fill=0 )
    if ( present(uBnd) ) then ! Store column sub-dimension bounds
      sparse%uBnd = uBnd
      if ( present(lBnd) ) then
        if ( size(ubnd,1) /= size(lbnd,1) ) &
          & call MLSMessage ( MLSMSG_Error, moduleName, "Size(LBnd) /= Size(UBnd)" )
        sparse%lBnd = lBnd
      else
        sparse%lBnd = uBnd ! Just to allocate it
        sparse%lBnd(:) = 1
      end if
    else ! Column doesn't have subdimensions
      sparse%uBnd = [nc]
      sparse%lBnd = [1]
    end if
    if ( present(ne) ) then
      allocate ( sparse%e(ne), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if
    if ( present(what) ) sparse%what = what

  end subroutine Create_Sparse

  subroutine Add_Element_Stru ( Sparse, E )
    ! Add an element with value E%V at row E%NR and column E&NC of Sparse
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    class(sparse_t), intent(inout) :: Sparse  ! The matrix
    type(sparse_element_t), intent(in) :: E  ! Element to add at (e%nr,e%nc)

    call add_element_value ( sparse, e%v, e%nr, e%nc )

  end subroutine Add_Element_Stru

  subroutine Add_Element_Value ( Sparse, V, R, C )
    ! Add an element with value V at row R and column C of Sparse
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate

    class(sparse_t), intent(inout) :: Sparse  ! The matrix
    real(rp), intent(in) :: V                ! The value of the element
    integer, intent(in) :: R                 ! The row of the element
    integer, intent(in) :: C                 ! The col of the element

    type(sparse_element_t), allocatable :: E(:)
    integer :: FC, FR, NC, NE, New, NR, Stat
    character(127) :: Msg

    if ( .not. allocated(sparse%e) ) then
      allocate ( sparse%e(2*size(sparse%rows)), stat=stat, errmsg=msg )
      call test_allocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if

    ne = size(sparse%e)
    if ( sparse%ne >= ne ) call resize_sparse ( sparse, 2*ne )

    sparse%ne = sparse%ne + 1
    new = sparse%ne
    sparse%e(new) = sparse_element_t ( v, r, c, 0, 0 )

    nr = sparse%rows(r)      ! Last column in row R
    if ( nr == 0 ) then
      fr = new
    else
      fr = sparse%e(nr)%nr   ! First item in row R is next after the last one
      sparse%e(nr)%nr = new  ! Next item after old last item is the new one
    end if
    sparse%e(new)%nr = fr    ! Next after new item in row R
    sparse%rows(r) = new     ! Last item in row R is the new item
    nc = sparse%cols(c)      ! Last row in column C
    if ( nc == 0 ) then
      fc = new
    else
      fc = sparse%e(nc)%nc   ! First row in column C is next after the last one
      sparse%e(nc)%nc = new  ! Next item after old last item is the new one
    end if
    sparse%e(new)%nc = fc    ! Next after new item in column C
    sparse%cols(c) = new     ! Last item in Column c is the new item

  end subroutine Add_Element_Value

  subroutine Add_Elements_Stru ( Sparse, E )
    ! Add elements with values E%V at rows E%NR and columns E%NC of Sparse
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    class(sparse_t), intent(inout) :: Sparse     ! The matrix
    type(sparse_element_t), intent(in) :: E(:)   ! Elements to add at (e%nr,e%nc)

    integer :: I

    do i = 1, size(e,1)
      call add_element_value ( sparse, e(i)%v, e(i)%nr, e(i)%nc )
    end do

  end subroutine Add_Elements_Stru

  subroutine Dump_Sparse ( Sparse, Name, Format, Colon, Row )
    use Array_Stuff, only: Subscripts
    use Dump_Options, only: SDFormatDefault
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String
    class(sparse_t), intent(in) :: Sparse         ! The matrix
    character(*), intent(in), optional :: Name    ! Print it if present
    character(*), intent(in), optional :: Format  ! Else sdFormatDefault
    logical, intent(in), optional :: Colon        ! Single subscript followed by
                                                  ! colon instead of within
                                                  ! parentheses, default false
    integer, intent(in), optional :: Row          ! Only dump this row

    integer :: CS(size(sparse%uBnd,1))
    integer :: I, J, K, L, R1, R2, S
    logical :: More  ! Can print more than the name
    logical :: MyColon
    character(len=64) :: MyFormat

    myFormat = sdFormatDefault
    if ( present(format) ) myFormat = format
    myColon = .false.
    if ( present(colon) ) myColon = colon

    ! Print the name and dimensions
    if ( present(name) ) call output ( name )
    if ( sparse%what /= 0 ) then
      if ( present(name) ) call blanks ( 1 )
      call display_string ( sparse%what )
    end if
    if ( present(name) .or. sparse%what /= 0 ) then
      call output ( size(sparse%rows), before='(' )
      if ( any(sparse%lBnd /= 1 .or. sparse%uBnd /= 1) ) then
        do s = 1, size(sparse%ubnd)
          if ( sparse%ubnd(s) - sparse%lbnd(s) > 0 ) then
            call output ( ',' )
            if ( sparse%lbnd(s) /= 1 ) call output ( sparse%lbnd(s), after=':' )
            call output ( sparse%ubnd(s) )
          end if
        end do
      else
        call output ( size(sparse%cols), before=',' )
      end if
      call output ( ')', advance='yes' )
    end if

    ! Print the nonzero array elements, their row numbers and their column
    ! sub-dimension subscripts
    more = .true. ! Assume sparse%rows and sparse%e are allocated
    if ( .not. allocated(sparse%rows) ) then
      call output ( 'Sparse%Rows not allocated', advance='yes' )
      more = .false.
    end if
    if ( .not. allocated(sparse%e) ) then
      call output ( 'Sparse%E not allocated', advance='yes' )
      more = .false.
    end if
    if ( more ) then
      r1 = 1
      r2 = size(sparse%rows,1)
      if ( present(row) ) then
        r1 = row
        r2 = row
      end if
      do i = r1, r2
        j = sparse%rows(i)                     ! First column
        if ( j == 0 ) cycle                    ! Don't dump an empty row
        call output ( i, places=4, after="#" ) ! Row number
        k = 0
        do
          if ( k == 5 ) then
            call newLine
            call blanks ( 5 )
          end if
          j = sparse%e(j)%nr                   ! Next column in this row
          k = k + 1
          cs = subscripts ( sparse%e(j)%c, sparse%uBnd, sparse%lBnd )
          if ( size(cs) == 1 .and. myColon ) then
            call output ( cs(1), places=4, after=": " )
          else
            call output ( cs(1), before=" (" )
            do l = 2, size(cs)
              call output ( cs(l), before="," )
            end do
            call output ( ") " )
          end if
          call output ( sparse%e(j)%v, format=myFormat )
          if ( j == sparse%rows(i) ) exit      ! No more columns this row
        end do ! Columns
        call newLine
      end do ! Rows
    end if

  end subroutine Dump_Sparse

  subroutine Destroy_Sparse ( Sparse )
    ! Create a sparse matrix representation
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    class(sparse_t), intent(inout) :: Sparse  ! The matrix

    character(127) :: Msg
    integer :: Stat

    sparse%ne = 0
    call deallocate_test ( sparse%rows, "Sparse%Rows", moduleName )
    call deallocate_test ( sparse%cols, "Sparse%Cols", moduleName )
    if ( allocated(sparse%e) ) then
      deallocate ( sparse%e, stat=stat, errmsg=msg )
      call test_deallocate ( stat, moduleName, "Sparse%E", ermsg=msg )
    end if

  end subroutine Destroy_Sparse

  pure real(rp) function Row_Dot_Vec ( Sparse, R, Vector ) result ( D )
    ! Compute dot product of row R of Sparse with Vector
    class(sparse_t), intent(in) :: Sparse  ! The matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(in) :: Vector(:)      ! The vector

    integer :: NC                          ! Next column in the row

    d = 0.0
    nc = sparse%rows(r)                    ! Last element in the row
    if ( nc /= 0 ) then                    ! Any columns in this row?
      do
        d = d + sparse%e(nc)%v * vector(sparse%e(nc)%c)
        nc = sparse%e(nc)%nc               ! Next column in the row
        if ( nc == sparse%rows(r) ) exit   ! back to the last one?
      end do
    end if

  end function Row_Dot_Vec

  subroutine Resize_Sparse ( Sparse, NE )
    ! Resize Sparse%E with new size NE
    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate, Test_Deallocate
    class(sparse_t), intent(inout) :: Sparse  ! The matrix
    integer, intent(in), optional :: NE       ! New size for Sparse%E,
                                              ! Default Sparse%NE

    type(sparse_element_t), allocatable :: E(:)
    integer :: N, Stat
    character(127) :: Msg

    n = sparse%ne
    if ( present(ne) ) n = ne
    allocate ( e(n), stat=stat, errmsg=msg )
    call test_allocate ( stat, moduleName, "E", ermsg=msg )
    e(1:sparse%ne) = sparse%e(1:sparse%ne)
    deallocate ( sparse%e, stat=stat, errmsg=msg )
    call test_deallocate ( stat, moduleName, "E", ermsg=msg )
    call move_alloc ( e, sparse%e )

  end subroutine Resize_Sparse

  subroutine Sparse_Dot_Vec ( Sparse, Vector, Prod )
    ! Multiply Sparse by Vector producing Product
    class(sparse_t), intent(in) :: Sparse  ! The matrix
    real(rp), intent(in) :: Vector(:)      ! The vector Sparse is to multiply
    real(rp), intent(out) :: Prod(:)

    integer :: R

    do r = 1, size(sparse%rows)
      prod(r) = row_dot_vec ( sparse, r, vector )
    end do

  end subroutine Sparse_Dot_Vec

  pure real(rp) function Vec_Dot_Col ( Vector, Sparse, C ) result ( D )
    ! Compute dot product of Vector with column C of Sparse
    real(rp), intent(in) :: Vector(:)      ! The vector
    class(sparse_t), intent(in) :: Sparse  ! The matrix
    integer, intent(in) :: C               ! Which col to multiply by Vector

    integer :: NR                          ! Next row in the column

    d = 0.0
    nr = sparse%cols(c)                    ! Last element in the column
    if ( nr /= 0 ) then                    ! Any rows in this column?
      do
        nr = sparse%e(nr)%nr               ! Next row in the column
        d = d + sparse%e(nr)%v * vector(sparse%e(nr)%r)
        if ( nr == sparse%cols(c) ) exit   ! back to the last one?
      end do
    end if

  end function Vec_Dot_Col

  subroutine Vec_Dot_Sparse ( Vector, Sparse, Prod )
    ! Multiply Sparse by Vector producing Product
    real(rp), intent(in) :: Vector(:)      ! The vector Sparse is to multiply
    class(sparse_t), intent(in) :: Sparse  ! The matrix
    real(rp), intent(out) :: Prod(:)

    integer :: C

    do c = 1, size(sparse%cols)
      prod(c) = vec_dot_col ( vector, sparse, c )
    end do

  end subroutine Vec_Dot_Sparse

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Sparse_m

! $Log$
! Revision 2.1  2017/11/01 18:52:11  vsnyder
! Initial commit
!
