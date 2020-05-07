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
  public :: Sparse_t, Sparse_Element_t

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
    integer :: NRows = 0            ! How many rows are filled
    integer, allocatable :: LBnd(:) ! Column sub-dimension lower bounds
    integer, allocatable :: UBnd(:) ! Column sub-dimension upper bounds
    type(sparse_element_t), allocatable :: E(:) ! nonzero elements
  contains
    procedure :: Add_Elements_Stru
    procedure :: Add_Element_Stru
    procedure :: Add_Element_Value
    generic :: Add_Element => Add_Element_Stru, Add_Element_Value, Add_Elements_Stru
    procedure :: Check_Column_Dup
    procedure :: Clear_All_Flags => Sparse_Clear_All_Flags
    procedure :: Clear_Col_Vec => Sparse_Clear_Col ! Clear a column vector
    procedure :: Clear_Col_And_Flags => Sparse_Clear_Col_And_Flags
    generic :: Clear_Col => Clear_Col_Vec, Clear_Col_And_Flags
    procedure :: Clear_Col_Flags => Sparse_Clear_Col_Flags
    procedure :: Clear_Vec => Sparse_Clear_Vec ! Clear a row vector
    generic :: Clear_Flags => Clear_All_Flags, Clear_Col_Flags
    procedure :: Copy => Copy_Sparse
    procedure :: Create => Create_Sparse
    procedure :: Densify
    procedure :: Destroy => Destroy_Sparse
    procedure :: Dump => Dump_Sparse
    procedure :: Empty => Empty_Sparse
    procedure :: Get_All_Flags => Sparse_Get_All_Flags
    procedure :: Get_Col_Flags => Sparse_Get_Col_Flags
    procedure :: Get_Col_Vec => Sparse_Get_Col
    procedure :: Get_Col_Vec_And_Flags => Sparse_Get_Col_And_Flags
    procedure :: Get_Col_Vec_And_Sparsity => Sparse_Get_Col_And_Sparsity
    generic :: Get_Col => Get_Col_Vec, Get_Col_Vec_And_Flags, &
             & Get_Col_Vec_And_Sparsity
    generic :: Get_Flags => Get_All_Flags, Get_Col_Flags
    procedure :: Invert_Column_Indices
    procedure :: List => List_Sparse
    procedure :: Resize => Resize_Sparse
    procedure :: Row_Dot_Matrix
    procedure :: Row_Dot_Vec_1D
    procedure :: Row_Dot_Vec_2D
    generic :: Row_Dot_Vec => Row_Dot_Vec_1D, Row_Dot_Vec_2D
    procedure :: Row_Times_Vec
    procedure :: Sparse_Dot_Vec_1D
    procedure :: Sparse_Dot_Vec_2D
    generic :: Sparse_Dot_Vec => Sparse_Dot_Vec_1D, Sparse_Dot_Vec_2D
    procedure, pass(sparse) :: Vec_Dot_Col
    procedure, pass(sparse) :: Vec_Dot_Sparse
  end type Sparse_t

  interface Add_Element
    module procedure Add_Element_Stru, Add_Element_Value, Add_Elements_Stru
  end interface Add_Element

  interface Allocate_Test
    module procedure Allocate_Test_Elements
  end interface

  interface DeAllocate_Test
    module procedure DeAllocate_Test_Elements
  end interface

  interface Dump
    module procedure Dump_Sparse, Dump_Sparse_Array
  end interface Dump

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------- -----------------

contains

  ! ----------------------------------------------  Create_Sparse  -----
  subroutine Create_Sparse ( Sparse, NR, NC, NE, UBnd, LBnd, What )
    ! Create a sparse matrix representation
    use Allocate_Deallocate, only: Allocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    class(sparse_t), intent(inout) :: Sparse  ! The sparse matrix
    integer, intent(in) :: NR                 ! Number of rows
    integer, intent(in) :: NC                 ! Number of cols
    integer, intent(in), optional :: NE       ! Initial size ( sparse%e )
    integer, intent(in), optional :: UBnd(:)  ! Column sub-dimensions upper bounds
    integer, intent(in), optional :: LBnd(:)  ! Column sub-dimensions lower bounds
                                              ! Not used if UBnd not present
    integer, intent(in), optional :: What     ! A string index, for dumps

    logical :: ReCreate ! Allocate or reallocate

    reCreate = .not. allocated(sparse%rows)
    if ( .not. reCreate ) reCreate = size(sparse%rows) /= nr
    if ( reCreate ) then
      call allocate_test ( sparse%rows, nr, "Sparse%Rows", moduleName, fill=0 )
    else
      sparse%rows = 0
    end if
    reCreate = .not. allocated(sparse%cols)
    if ( .not. reCreate ) reCreate = size(sparse%cols) /= nc
    if ( reCreate ) then
      call allocate_test ( sparse%cols, nc, "Sparse%Cols", moduleName, fill=0 )
    else
      sparse%cols = 0
    end if
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
    else ! Column doesn't have sub-dimensions
      sparse%uBnd = [nc]
      sparse%lBnd = [1]
    end if
    if ( present(ne) ) then
      reCreate = .not. allocated(sparse%e)
      if ( .not. reCreate ) reCreate = size(sparse%e) < ne
      if ( reCreate ) then
        call deallocate_test ( sparse%e, "Sparse%E in Create_Sparse" )
        call allocate_test ( sparse%e, ne, "Sparse%E in Create_Sparse" )
      end if
    end if
    sparse%ne = 0
    sparse%nRows = 0
    if ( present(what) ) sparse%what = what

  end subroutine Create_Sparse

  ! -------------------------------------------  Add_Element_Stru  -----
  subroutine Add_Element_Stru ( Sparse, E )
    ! Add an element with value E%V at row E%NR and column E&NC of Sparse
    class(sparse_t), intent(inout) :: Sparse  ! The sparse matrix
    type(sparse_element_t), intent(in) :: E  ! Element to add at (e%nr,e%nc)

    call add_element_value ( sparse, e%v, e%nr, e%nc )

  end subroutine Add_Element_Stru

  ! ------------------------------------------  Add_Element_Value  -----
  subroutine Add_Element_Value ( Sparse, V, R, C )
    ! Add an element with value V at row R and column C of Sparse
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    class(sparse_t), intent(inout) :: Sparse ! The sparse matrix
    real(rp), intent(in) :: V                ! The value of the element
    integer, intent(in) :: R                 ! The row of the element
    integer, intent(in) :: C                 ! The col of the element

    integer :: FC, FR, NC, NE, New, NR

    if ( v == 0 ) return ! Don't add zero elements!

    if ( .not. allocated(sparse%e) ) then
      if ( .not. allocated(sparse%rows) ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, "In Add_Element_Value, Sparse not yet created" )
      call allocate_test ( sparse%e, 2*size(sparse%rows), &
        & "Sparse%E in Add_Element_Value" )
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

  ! ------------------------------------------  Add_Elements_Stru  -----
  subroutine Add_Elements_Stru ( Sparse, E )
    ! Add elements with values E%V at rows E%NR and columns E%NC of Sparse
    class(sparse_t), intent(inout) :: Sparse     ! The sparse matrix
    type(sparse_element_t), intent(in) :: E(:)   ! Elements to add at (e%nr,e%nc)

    integer :: I

    do i = 1, size(e,1)
      call add_element_value ( sparse, e(i)%v, e(i)%nr, e(i)%nc )
    end do

  end subroutine Add_Elements_Stru

  ! -------------------------------------  Allocate_Test_Elements  -----
  subroutine Allocate_Test_Elements ( Elements, Size, What )
    use Allocate_Deallocate, only: Test_Allocate
    type(sparse_element_t), allocatable, intent(inout) :: Elements(:)
    integer, intent(in) :: Size
    character(*), intent(in) :: What
    integer :: Stat
    character(127) :: Msg
    call deallocate_test ( elements, what )
    allocate ( elements(size), stat=stat, errmsg=msg )
    call test_allocate ( stat, moduleName, what, ermsg=msg )
  end subroutine Allocate_Test_Elements

  ! -------------------------------------------  Check_Column_Dup  -----
  subroutine Check_Column_Dup ( Sparse, Name, Found )
    use Output_m, only: Output
    class(sparse_t), intent(in) :: Sparse         ! The sparse matrix
    character(*), intent(in) :: Name
    logical, intent(out), optional :: Found       ! Found a duplicate
    integer :: C, R0, R1, R2
    logical :: Named
    named = .false.
    if ( present(found) ) found=.false.
    do c = 1, size(sparse%cols)
      r0 = sparse%cols(c)
      if ( r0 == 0 ) cycle
      r1 = r0
      do
        r1 = sparse%e(r1)%nc ! Next element in same column
        if ( r1 == r0 ) exit
        r2 = sparse%e(r1)%nc
        do
          r2 = sparse%e(r2)%nc
          if ( r1 /= r2 ) then
            if ( sparse%e(r1)%r == sparse%e(r2)%r .and. &
               & sparse%e(r1)%c == sparse%e(r2)%c ) then
              if ( .not. named ) call output ( 'In sparse matrix ' // trim(name), &
                & advance='yes' )
              named=.true.
              if ( present(found) ) found=.true.
              call output ( c, before='In column ' )
              call output ( r1, before=', elements ' )
              call output ( r2, before=' and ' )
              call output ( sparse%e(r1)%r, before=' are duplicates for row ' )
              if  ( sparse%e(r1)%v == sparse%e(r2)%v ) then
                call output ( '.  Values are the same.' , advance='yes' )
              else
                call output ( '.  Values differ.' , advance='yes' )
              end if
            end if
          end if
          if ( r2 == r0 ) exit
        end do
      end do
    end do
  end subroutine Check_Column_Dup

  ! ------------------------------------------------  Copy_Sparse  -----
  subroutine Copy_Sparse ( To_Sparse, From_Sparse )
  ! Copy From_Sparse to To_Sparse without fiddling with allocations
  ! any more than necessary
    use Allocate_Deallocate, only: Deallocate_Test
    class(sparse_t), intent(inout) :: To_Sparse  ! The "To" sparse matrix
    class(sparse_t), intent(in) :: From_Sparse  ! The "From" sparse matrix

    if ( allocated(from_sparse%rows) ) then
      to_sparse%rows = from_sparse%rows
    else
      call deallocate_test ( to_sparse%rows, 'To_Sparse%rows', moduleName )
    end if

    if ( allocated(from_sparse%cols) ) then
      to_sparse%cols = from_sparse%cols
    else
      call deallocate_test ( to_sparse%cols, 'To_Sparse%Cols', moduleName )
    end if

    if ( allocated(from_sparse%lbnd) ) then
      to_sparse%lbnd = from_sparse%lbnd
    else
      call deallocate_test ( to_sparse%lbnd, 'To_Sparse%Lbnd', moduleName )
    end if

    if ( allocated(from_sparse%ubnd) ) then
      to_sparse%ubnd = from_sparse%ubnd
    else
      call deallocate_test ( to_sparse%ubnd, 'To_Sparse%Ubnd', moduleName )
    end if

    to_sparse%ne = from_sparse%ne
    to_sparse%nRows = from_sparse%nRows
    to_sparse%what = from_sparse%what

    if ( allocated(from_sparse%e) ) then
      if ( allocated(to_sparse%e) ) then
        if ( size(to_sparse%e) < size(from_sparse%e) ) then
          to_sparse%e = from_sparse%e ! Re-allocates to_sparse%e
        else
          to_sparse%e(:size(from_sparse%e)) = from_sparse%e
        end if
      else
        to_sparse%e = from_sparse%e
      end if
    else
      call deallocate_test ( to_sparse%e, 'To_Sparse%E in Copy_Sparse' )
    end if

  end subroutine Copy_Sparse

  ! -----------------------------------  DeAllocate_Test_Elements  -----
  subroutine DeAllocate_Test_Elements ( Elements, What )
    use Allocate_Deallocate, only: Test_DeAllocate
    type(sparse_element_t), allocatable, intent(inout) :: Elements(:)
    character(*), intent(in) :: What
    integer :: Stat
    character(127) :: Msg
    if ( allocated(elements) ) then
      deallocate ( elements, stat=stat, errmsg=msg )
      call test_deallocate ( stat, moduleName, what, ermsg=msg )
    end if
  end subroutine DeAllocate_Test_Elements

  ! ----------------------------------------------------  Densify  -----
  subroutine Densify ( Sparse, Dense )
    ! Copy Sparse to a dense representation in Dense, probably for dumping
    use Allocate_Deallocate, only: Allocate_test
    class(sparse_t), intent(in) :: Sparse            ! The sparse matrix
    real(rp), allocatable, intent(out) :: Dense(:,:) ! Dense representation
    integer :: J ! Column subscript
    if ( .not. allocated(sparse%rows) .or. .not. allocated(sparse%cols) ) return
    call allocate_test ( dense, size(sparse%rows,1), size(sparse%cols,1), &
      & 'Dense', moduleName, fill=0.0_rp )
    do j = 1, size(sparse%cols)
      call sparse%get_col_vec ( j, dense(:,j) )
    end do
  end subroutine Densify

  ! ------------------------------------------------  Dump_Sparse  -----
  subroutine Dump_Sparse ( Sparse, Name, Format, Colon, This, Offset, Transpose, &
                         & Width, Exclude, One_D_Col )
    use Array_Stuff, only: Subscripts
    use Dump_Options, only: SDFormatDefault
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String
    class(sparse_t), intent(in) :: Sparse         ! The sparse matrix
    character(*), intent(in), optional :: Name    ! Print it if present
    character(*), intent(in), optional :: Format  ! Else sdFormatDefault
    logical, intent(in), optional :: Colon        ! Single subscript followed by
                                                  ! colon instead of within
                                                  ! parentheses, default false
    integer, intent(in), optional :: This         ! Only dump this row or column
                                                  ! (depending on transpose)
    logical, intent(in), optional :: Offset       ! Print column numbers starting
                                                  ! at 1 instead of their positions
                                                  ! in Grids_t
    logical, intent(in), optional :: Transpose    ! Dump transpose -- each row
                                                  ! of output is a column
    integer, intent(in), optional :: Width        ! How many per line, default 5
    real(rp), intent(in), optional :: Exclude     ! Don't print these
    logical, intent(in), optional :: One_D_Col    ! Print 1D column subscript

    integer :: CS(size(sparse%uBnd,1))
    integer :: C1, C2, I, J, K, L, R1, R2, S
    logical :: DidOne
    logical :: DoThis ! Do this row or column
    logical :: More   ! Can print more than the name
    logical :: MyColon
    character(len=64) :: MyFormat
    logical :: MyOffset
    logical :: MyOne_D_Col
    logical :: MyTranspose
    integer :: MyWidth
    integer :: Places

    if ( .not. allocated(sparse%rows) .or. .not. allocated(sparse%cols) ) &
      call MLSMessage ( MLSMSG_Error, moduleName, "Can't dump Eta that's not created" )

    myColon = .false.
    if ( present(colon) ) myColon = colon
    myFormat = sdFormatDefault
    if ( present(format) ) myFormat = format
    MyOffset = .false.
    if ( present(offset) ) MyOffset = offset
    myOne_d_col = .false.
    if ( present(one_d_col) ) myOne_d_col = one_d_col
    MyTranspose = .false.
    if ( present(Transpose) ) MyTranspose = transpose !.and. size(cs,1) > 1
    myWidth = 5
    if ( present(width) ) myWidth = width

    ! Print the name and dimensions
    if ( present(name) ) call output ( name )
    if ( sparse%what /= 0 ) then
      if ( present(name) ) call blanks ( 1 )
      call display_string ( sparse%what )
    end if
    if ( present(name) .or. sparse%what /= 0 ) then
      call output ( sparse%nRows, before='(' )
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

    more = .true. ! Assume sparse%rows and sparse%e are allocated
    if ( .not. allocated(sparse%rows) ) then
      call output ( 'Sparse%Rows not allocated' )
      more = .false.
    end if
    if ( .not. allocated(sparse%e) ) then
      call output ( 'Sparse%E not allocated' )
      more = .false.
    end if
    if ( sparse%ne <= 0 ) then
      call output ( ' Sparse%NE <= 0' )
      more=.false.
    end if
    didOne = .false.
    if ( more ) then
      if ( myTranspose ) then
        ! Print the nonzero array elements, their column sub-dimension
        ! subscripts and their row numbers.
        c1 = 1
        c2 = size(sparse%cols,1)
        if ( present(this) ) then
          c1 = this
          c2 = this
        end if
        do i = c1, c2
          j = sparse%cols(i)                     ! Last row this column
          if ( j == 0 ) cycle                    ! Don't dump an empty column
          doThis = .not. present(exclude)
          if ( .not. doThis ) then
            do
              j = sparse%e(j)%nc                 ! Next row in this column
              doThis = sparse%e(j)%v /= exclude
              if ( doThis ) exit
              if ( j == sparse%cols(i) ) exit    ! No more columns this row
            end do
            j = sparse%cols(i)
          end if
          if ( .not. doThis ) cycle
          didOne = .true.
          if ( myOne_d_col ) then
            call output ( sparse%e(j)%c, before=' (', after=') ' )
          else
            cs = subscripts ( sparse%e(j)%c, sparse%uBnd, sparse%lBnd )
            if ( myOffset ) cs = cs - sparse%lBnd + 1
            call output ( cs(1), before=" (" )
            do l = 2, size(cs)
              call output ( cs(l), before="," )
            end do
            if ( myColon ) then
              call output ( "):" )
            else
              call output ( ") " )
            end if!ma
          end if
          places = sum(int(log10(real(cs)))) + size(cs) + 3
          k = 0
          do
            if ( k >= myWidth ) then
              call newLine
              call blanks ( places )
              k = 0
            end if
            j = sparse%e(j)%nc                   ! Next row in this column
            k = k + 1
            call output ( sparse%e(j)%r, places=4 )
            if ( myColon ) call output ( ':' )
            call output ( sparse%e(j)%v, format=myFormat )
            if ( j == sparse%cols(i) ) exit      ! No more columns this row
          end do ! Rows
          call newLine
        end do ! Columns
      else ! .not. myTranspose
        ! Print the nonzero array elements, their row numbers and their
        ! column sub-dimension subscripts.
        r1 = 1
        r2 = sparse%nRows
        if ( r2 == 0 ) r2 = size(sparse%rows,1)
        if ( present(this) ) then
          r1 = this
          r2 = this
        end if
        do i = r1, r2
          j = sparse%rows(i)                     ! Last column this row
          if ( j == 0 ) cycle                    ! Don't dump an empty row
          doThis = .not. present(exclude)
          if ( .not. doThis ) then
            do
              j = sparse%e(j)%nr                 ! Next column in this row
              doThis = sparse%e(j)%v /= exclude
              if ( doThis ) exit
              if ( j == sparse%rows(i) ) exit    ! No more rows this column
            end do
            j = sparse%rows(i)
          end if
          if ( .not. doThis ) cycle
          didOne = .true.
          if ( sparse%nRows > 1 ) &
            & call output ( i, places=4, after="#" ) ! Row number
          k = 0
          do
            if ( k >= myWidth ) then
              call newLine
              call blanks ( 5 )
              k = 0
            end if
            j = sparse%e(j)%nr                   ! Next column in this row
            k = k + 1
            if ( myOne_d_col ) then
              call output ( sparse%e(j)%c, before=' (', after=') ' )
            else
              cs = subscripts ( sparse%e(j)%c, sparse%uBnd, sparse%lBnd )
              if ( myOffset ) cs = cs - sparse%lBnd + 1
              if ( size(cs) == 1 .and. myColon ) then
                call output ( cs(1), places=4, after=": " )
              else
                call output ( cs(1), before=" (" )
                do l = 2, size(cs)
                  call output ( cs(l), before="," )
                end do
                if ( myColon ) then
                  call output ( "): " )
                else
                  call output ( ") " )
                end if
              end if
            end if
            call output ( sparse%e(j)%v, format=myFormat )
            if ( j == sparse%rows(i) ) exit      ! No more columns this row
          end do ! Columns
          call newLine
        end do ! Rows
      end if
    else
      call newLine
    end if

    if ( .not. didOne ) then
      if ( present(exclude) ) then
        call output ( exclude, &
          & before="Either there are no elements, or they're all = ", &
          & advance="yes" )
      else
        call output ( "There are no nonzero elements", advance="yes" )
      end if
    end if

  end subroutine Dump_Sparse
  
  ! ------------------------------------------  Dump_Sparse_Array  -----
  subroutine Dump_Sparse_Array ( Sparse, Name, Format, Colon, This, Offset, &
                               & Transpose, Width )
    class(sparse_t), intent(in) :: Sparse(:)      ! The sparse matrix
    character(*), intent(in), optional :: Name    ! Print it if present
    character(*), intent(in), optional :: Format  ! Else sdFormatDefault
    logical, intent(in), optional :: Colon        ! Single subscript followed by
                                                  ! colon instead of within
                                                  ! parentheses, default false
    integer, intent(in), optional :: This         ! Only dump this row or column
                                                  ! (depending on transpose)
    logical, intent(in), optional :: Offset       ! Print column numbers starting
                                                  ! at 1 instead of their positions
                                                  ! in Grids_t
    logical, intent(in), optional :: Transpose    ! Dump transpose -- each row
                                                  ! of output is a column
    integer, intent(in), optional :: Width        ! How many per line, default 5

    integer :: Sps

    do sps = 1, size(sparse)
      call sparse(sps)%dump ( name, format, colon, this, offset, transpose, width )
    end do
  end subroutine Dump_Sparse_Array

  ! ---------------------------------------------  Destroy_Sparse  -----
  subroutine Destroy_Sparse ( Sparse_Not_Used )
    ! Destroy a sparse matrix representation by default initializing it,
    ! which deallocates any of its allocated allocatable components.
    class(sparse_t), intent(out) :: Sparse_Not_Used  ! The sparse matrix
  end subroutine Destroy_Sparse

  ! -----------------------------------------------  Empty_Sparse  -----
  subroutine Empty_Sparse ( Sparse, Resize )
    ! Set Sparse%ne, sparse%rows, and sparse%cols to zero to make it effectively
    ! empty, but without deallocating everything.  Useful if the next time it's
    ! used it will have the same number (or fewer) of rows.  If the number of
    ! columns might be different, it ought to be destroyed because so many
    ! things depend upon the size of Sparse%Cols.
    class(sparse_t), intent(inout) :: Sparse
    logical, intent(in), optional :: Resize ! Resize Sparse%E to
                                            ! 2*size(sparse%rows)
    sparse%ne = 0
    if ( allocated(sparse%rows) ) sparse%rows = 0
    if ( allocated(sparse%cols) ) sparse%cols = 0
    if ( present(resize) ) then
      if ( resize ) call allocate_test ( sparse%e, 2*size(sparse%rows), &
        & "Sparse%E in Empty_Sparse" )
    end if
  end subroutine Empty_Sparse

  ! --------------------------------------  Invert_Column_Indices  -----
  subroutine Invert_Column_Indices ( Sparse )
    ! Invert the column indices because the column basis was inverted.
    ! The routine to compute Eta only works with increasing basis, so the
    ! basis had to be reversed on the call, which made the column indices
    ! inverted.
    class(sparse_t), intent(inout) :: Sparse
    integer :: R, C, L
    do r = 1, sparse%nRows
      l = sparse%rows(r) ! Last column in row R
      c = sparse%e(l)%nr ! First element in row R
      do
        sparse%e(c)%c = size(sparse%cols) - sparse%e(c)%c + 1
        if ( c == l ) exit
        c = sparse%e(c)%nr
      end do
    end do
  end subroutine Invert_Column_Indices

  ! ------------------------------------------------  List_Sparse  -----
  subroutine List_Sparse ( Sparse, Name, Details )
    ! List the contents of Sparse_t without trying to follow the links for
    ! rows and columns.  Useful for debugging if the structure is messed up.
    use Dump_0, only: Dump
    use Output_m, only: NewLine, Output
    use String_Table, only: Display_String
    class(sparse_t), intent(in) :: Sparse
    character(*), intent(in), optional :: Name
    integer, intent(in), optional :: Details ! <= 0 => Only the structure
                                             ! > 0 => The elements too.
                                             ! Default 1
    integer :: F, L ! Locations of first, last nonzeros in sparse%col
    integer :: I, MyDetails

    if ( present(name) ) call output ( name )
    if ( sparse%what /= 0 ) then
      if ( present(name) ) call output ( ' ' )
      call display_string ( sparse%what )
    end if
    if ( present(name) .or. sparse%what /= 0 ) call newLine
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( .not. allocated(sparse%rows) ) then
      call output ( 'Sparse matrix not initialized', advance='yes' )
      return
    end if
    call dump ( sparse%rows(1:sparse%nRows), name='Rows' )
    do f = 1, size(sparse%cols)
      if ( sparse%cols(f) /= 0 ) exit
    end do
    do l = size(sparse%cols), 1, -1
      if ( sparse%cols(l) /= 0 ) exit
    end do
    call dump ( sparse%cols(f:l), name='Cols', lbound=f )
    call dump ( sparse%lbnd, name='Lbnd' )
    call dump ( sparse%ubnd, name='Ubnd' )
    call output ( sparse%nRows, before='NRows ', advance='yes' )
    if ( sparse%ne /= 0 .and. myDetails > 0 ) then
      call output ( sparse%ne, before='         R    C   NR   NC   V  NE = ', &
                  & advance='yes' )
      do i = 1, sparse%ne
        call output ( i, places=4, after='#' )
        call output ( sparse%e(i)%r, places=5 )
        call output ( sparse%e(i)%c, places=5 )
        call output ( sparse%e(i)%nr, places=5 )
        call output ( sparse%e(i)%nc, places=5 )
        call output ( sparse%e(i)%v, before='  ', advance='yes' )
      end do
    else if ( sparse%ne == 0 ) then
      call output ( 'NE is zero', advance='yes' )
    end if
  end subroutine List_Sparse

  ! ----------------------------------------------  Resize_Sparse  -----
  subroutine Resize_Sparse ( Sparse, NE )
    ! Resize Sparse%E with new size NE
    class(sparse_t), intent(inout) :: Sparse  ! The sparse matrix
    integer, intent(in), optional :: NE       ! New size for Sparse%E,
                                              ! Default Sparse%NE

    type(sparse_element_t), allocatable :: E(:)
    integer :: N

    n = sparse%ne
    if ( present(ne) ) n = max(ne,n)
    call allocate_test ( e, n, "E in Resize_Sparse" )
    e(1:sparse%ne) = sparse%e(1:sparse%ne)
    call deallocate_test ( sparse%e, "Sparse%E in Resize_Sparse" )
    call move_alloc ( e, sparse%e )

  end subroutine Resize_Sparse

  ! ---------------------------------------------  Row_Dot_Matrix  -----
  pure subroutine Row_Dot_Matrix ( Sparse, R, Matrix, Prod, ExpLog, LogOnly, &
                                 & Rows )
    ! Compute vector-matrix product, where the vector is row R of Sparse.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(in) :: Matrix(:,:)    ! The matrix
    real(rp), intent(out) :: Prod(:)       ! Same extent as size(matrix,1)
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority
    class(sparse_t), intent(in), optional :: Rows ! It's only necessary to
                                           ! compute elements of the product
                                           ! where Rows%cols is not zero,
                                           ! usually because Prod will later
                                           ! be used in rows%row_dot_vec.

    integer :: J                           ! Column index in the row

    if ( present(rows) ) then
      do j = 1, size(matrix,1)
        if ( rows%cols(j) /= 0 ) &
          prod(j) = sparse%row_dot_vec ( r, matrix(j,:), expLog, logOnly )
      end do
    else
      do j = 1, size(matrix,1)
        prod(j) = sparse%row_dot_vec ( r, matrix(j,:), expLog, logOnly )
      end do
    end if

  end subroutine Row_Dot_Matrix

  ! ---------------------------------------------  Row_Dot_Vec_1D  -----
  pure real(rp) function Row_Dot_Vec_1D ( Sparse, R, Vector, ExpLog, LogOnly ) &
    & result ( D )
    ! Compute dot product of row R of Sparse with Vector
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(in) :: Vector(:)      ! The vector
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority

    integer :: J                           ! Column index in the row
    logical :: MyExpLog, MyLog

    myLog = .false.
    if ( present(logOnly) ) myLog = logOnly
    myExpLog = .false.
    if ( present(expLog) ) myExpLog = expLog
    myLog = myLog .or. myExpLog

    d = 0.0
    j = sparse%rows(r)                     ! Last element in the row
    if ( j /= 0 ) then                     ! Any columns in this row?
      if ( myLog ) then
        do
          d = d + sparse%e(j)%v * log(max(vector(sparse%e(j)%c),1.0e-9_rp))
          j = sparse%e(j)%nr                ! Next column in the row
          if ( j == sparse%rows(r) ) exit   ! back to the last one?
        end do
      else
        do
          d = d + sparse%e(j)%v * vector(sparse%e(j)%c)
          j = sparse%e(j)%nr                ! Next column in the row
          if ( j == sparse%rows(r) ) exit   ! back to the last one?
        end do
      end if
    end if
    if ( myExpLog ) d = exp(d)

  end function Row_Dot_Vec_1D

  ! ---------------------------------------------  Row_Dot_Vec_2D  -----
!  pure &
  real(rp) function Row_Dot_Vec_2D ( Sparse, R, Vector, ExpLog, LogOnly ) &
    & result ( D )
    ! Compute dot product of row R of Sparse with Vector
    use Array_Stuff, only: Subscripts
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(in) :: Vector(:,:)    ! The vector.  Sparse%e(j)%c is
                                           ! the array-element order index.
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority

    integer :: C                           ! Array-element order subscript for
                                           ! Vector
    integer :: J                           ! Column index in the row
    logical :: MyExpLog, MyLog
    integer :: S(2)                        ! Subscripts of Vector

    myLog = .false.
    if ( present(logOnly) ) myLog = logOnly
    myExpLog = .false.
    if ( present(expLog) ) myExpLog = expLog
    myLog = myLog .or. myExpLog

    d = 0.0
    j = sparse%rows(r)                     ! Last element in the row
    if ( j /= 0 ) then                     ! Any columns in this row?
      if ( myLog ) then
        do
          c = sparse%e(j)%c
          s = subscripts ( c, shape(vector) )
          d = d + sparse%e(j)%v * log(max(vector(s(1),s(2)),1.0e-9_rp))
          j = sparse%e(j)%nr                ! Next column in the row
          if ( j == sparse%rows(r) ) exit   ! back to the last one?
        end do
      else
        do
          c = sparse%e(j)%c
          s = subscripts ( c, shape(vector) )
          d = d + sparse%e(j)%v * vector(s(1),s(2))
          j = sparse%e(j)%nr                ! Next column in the row
          if ( j == sparse%rows(r) ) exit   ! back to the last one?
        end do
      end if
    end if
    if ( myExpLog ) d = exp(d)


  end function Row_Dot_Vec_2D

  ! ----------------------------------------------  Row_Times_Vec  -----
  subroutine Row_Times_Vec ( Sparse, R, Vector, Product )
    ! Store the product of the nonzero elements of row R of Sparse and
    ! corresponding elements of Vector in corresponding elements of Product.
    ! This is an element-by-element product.  Product is not initially made
    ! zero.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(in) :: Vector(:)      ! The vector
    real(rp), intent(inout) :: Product(:)  ! Elements corresponding to nonzeroes
                                           ! in row R are replaced
    integer :: J
    j = sparse%rows(r)  ! Last element in row R
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nr   ! Element in next column of row R
      product(sparse%e(j)%c) = sparse%e(j)%v * vector(sparse%e(j)%c)
      if ( j == sparse%rows(r) ) exit
    end do
  end subroutine Row_Times_Vec

  ! -------------------------------------  Sparse_Clear_All_Flags  -----
  subroutine Sparse_Clear_All_Flags ( Sparse, Flags )
    ! Make values of Flags that correspond to nonzeroes of Sparse false.
    class(sparse_t), intent(in) :: Sparse
    logical, intent(inout) :: Flags(:,:)

    integer :: E          ! Index of Sparse(sps_i)%e

    do e = 1, sparse%ne
      flags(sparse%e(e)%r,sparse%e(e)%c) = .false.
    end do
  end subroutine Sparse_Clear_All_Flags

  ! -------------------------------------------  Sparse_Clear_Col  -----
  subroutine Sparse_Clear_Col ( Sparse, C, Vector )
    ! Clear elements of Vector that correspond to nonzero elements of column C
    ! of Sparse.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    real(rp), intent(inout) :: Vector(:)   ! The vector
    integer :: J
    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      vector(sparse%e(j)%r) = 0
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Clear_Col

  ! ---------------------------------  Sparse_Clear_Col_And_Flags  -----
  subroutine Sparse_Clear_Col_And_Flags ( Sparse, C, Vector, Flags, Last )
    ! Clear elements of Vector that correspond to nonzero elements of column C
    ! of Sparse.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    real(rp), intent(inout) :: Vector(:)   ! The vector
    logical, intent(inout) :: Flags(:)     ! Set false where vector is set to
                                           ! zero here
    integer, intent(in), optional :: Last  ! Do not process rows > Last
    integer :: J, MyLast, R

    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    myLast = sparse%nRows
    if ( present(last) ) myLast = last
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      r = sparse%e(j)%r
      if ( r <= myLast ) then
        vector(r) = 0
        flags(r) = .false.
      end if
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Clear_Col_And_Flags

  ! -------------------------------------  Sparse_Clear_Col_Flags  -----
  subroutine Sparse_Clear_Col_Flags ( Sparse, C, Flags )
    ! Set elements of Flags that correspond to nonzero elements of column C
    ! of Sparse to be false.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    logical, intent(inout) :: Flags(:)     ! The flags
    integer :: J
    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      flags(sparse%e(j)%r) = .false.
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Clear_Col_Flags

  ! -------------------------------------------  Sparse_Clear_Vec  -----
!  pure &
  subroutine Sparse_Clear_Vec ( Sparse, R, Vector )
    ! Clear elements of Vector that correspond to nonzero elements of row R
    ! of Sparse.  This is used to clean up after Row_Times_Vec.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: R               ! Which row to multiply by Vector
    real(rp), intent(inout) :: Vector(:)   ! The vector
    integer :: J
    j = sparse%rows(r)     ! Last element in row R
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nr   ! Element in next column of row R
      vector(sparse%e(j)%c) = 0
      if ( j == sparse%rows(r) ) exit
    end do
  end subroutine Sparse_Clear_Vec

  ! ------------------------------------------  Sparse_Dot_Matrix  -----
!  pure &
  subroutine Sparse_Dot_Matrix ( Sparse, Matrix, Prod, ExpLog, LogOnly )
    ! Multiply Sparse by Matrix producing Prod
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    real(rp), intent(in) :: Matrix(:,:)    ! The matrix Sparse is to multiply
    real(rp), intent(out) :: Prod(:,:)
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority

    integer :: C, R

    do r = 1, size(sparse%rows)
      do c = 1, size(matrix,2)
        prod(r,c) = sparse%row_dot_vec ( r, matrix(:,c), expLog, logOnly )
      end do
    end do

  end subroutine Sparse_Dot_Matrix

  ! ------------------------------------------  Sparse_Dot_Vec_1D  -----
!  pure &
  subroutine Sparse_Dot_Vec_1D ( Sparse, Vector, Prod, ExpLog, LogOnly )
    ! Multiply Sparse by Vector producing Prod
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    real(rp), intent(in) :: Vector(:)      ! The vector Sparse is to multiply
    real(rp), intent(out) :: Prod(:)
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority

    integer :: I
    logical :: MyExpLog, MyLog
    integer :: R ! Row subscripts of Sparse%E and Prod

    myLog = .false.
    if ( present(logOnly) ) myLog = logOnly
    myExpLog = .false.
    if ( present(expLog) ) myExpLog = expLog
    myLog = myLog .or. myExpLog

    prod = 0
    if ( myLog ) then
      do i = 1, sparse%ne
        r = sparse%e(i)%r
        prod(r) = prod(r) + sparse%e(i)%v * log(max(vector(sparse%e(i)%c),1.0e-9_rp))
      end do
    else
      do i = 1, sparse%ne
        r = sparse%e(i)%r
        prod(r) = prod(r) + sparse%e(i)%v * vector(sparse%e(i)%c)
      end do
    end if
    if ( myExpLog ) prod = exp(prod)

  end subroutine Sparse_Dot_Vec_1D

  ! ------------------------------------------  Sparse_Dot_Vec_2D  -----
  subroutine Sparse_Dot_Vec_2D ( Sparse, Vector, Prod, ExpLog, LogOnly )
    ! Multiply Sparse by Vector producing Prod
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    real(rp), intent(in), target :: Vector(:,:) ! The vector Sparse is to
                                           ! multiply. Sparse%e(j)%c is the
                                           ! array-element order index.
    real(rp), intent(out), target :: Prod(:) ! The size of Prod should be same
                                           ! as the row dimension of Sparse
    logical, intent(in), optional :: ExpLog ! compute exp(dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp)))
    logical, intent(in), optional :: LogOnly ! compute dot product(sparse, 
                                           ! log(max(vector,1.0e-9_rp));
                                           ! ExpLog takes priority

    integer :: N                           ! Number of rows to use
    integer :: R
                                           ! to R

    n = sparse%nRows
    if ( n == 0 ) n = size(sparse%rows)
    do r = 1, min(n,ubound(prod,1))
      prod(r) = row_dot_vec_2d ( sparse, r, vector, expLog, logOnly )
    end do

  end subroutine Sparse_Dot_Vec_2D

  ! ---------------------------------------  Sparse_Get_All_Flags  -----
!  pure &
  subroutine Sparse_Get_All_Flags ( Sparse, Flags )
    ! Make values of Flags that correspond to nonzeroes of Sparse true.
    class(sparse_t), intent(in) :: Sparse
    logical, intent(inout) :: Flags(:,:)

    integer :: E          ! Index of Sparse(sps_i)%e

    do e = 1, sparse%ne
      flags(sparse%e(e)%r,sparse%e(e)%c) = .true.
    end do

  end subroutine Sparse_Get_All_Flags

  ! ---------------------------------------------  Sparse_Get_Col  -----
!  pure &
  subroutine Sparse_Get_Col ( Sparse, C, Vector )
    ! Get elements of Vector that correspond to nonzero elements of column C
    ! of Sparse.  Vector is not initially made zero.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which row to multiply by Vector
    real(rp), intent(inout) :: Vector(:)   ! The vector
    integer :: J
    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nc   ! Element in next row of column c
      vector(sparse%e(j)%r) = sparse%e(j)%v
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Get_Col

  ! -----------------------------------  Sparse_Get_Col_And_Flags  -----
!  pure &
  subroutine Sparse_Get_Col_And_Flags ( Sparse, C, Vector, Flags, Last )
    ! Get elements of Vector that correspond to nonzero elements of column C
    ! of Sparse.  Vector is not initially made zero.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    real(rp), intent(inout) :: Vector(:)   ! The vector
    logical, intent(inout) :: Flags(:)     ! True where Vector gets a value;
                                           ! Not initially set false here
    integer, intent(in), optional :: Last  ! Do not get rows beyond Last
    integer :: J, MyLast, R

    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    myLast = sparse%nRows
    if ( present(last) ) myLast = last
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      r = sparse%e(j)%r
      if ( r <= myLast ) then
        vector(r) = sparse%e(j)%v
        flags(r) = .true.
      end if
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Get_Col_And_Flags

  ! --------------------------------  Sparse_Get_Col_And_Sparsity  -----
!  pure &
  subroutine Sparse_Get_Col_And_Sparsity ( Sparse, C, Vector, NNZ, NZ )
    ! Get elements of Vector that correspond to nonzero elements of column C
    ! of Sparse.  Vector is not initially made zero.  Set NNZ to the number of
    ! nonzeroes.  Set NZ to the row numbers of nonzeroes
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    real(rp), intent(inout) :: Vector(:)   ! The vector
    integer, intent(out) :: NNZ            ! Number of nonzeroes
    integer, intent(out) :: NZ(:)
    integer :: J
    j = sparse%cols(c)     ! Last element in column c
    nnz = 0
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      vector(sparse%e(j)%r) = sparse%e(j)%v
      nnz = nnz + 1
      nz(nnz) = sparse%e(j)%r
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Get_Col_And_Sparsity

  ! -------------------------------------  Sparse_Get_Col_Flags  -----
!  pure &
  subroutine Sparse_Get_Col_Flags ( Sparse, C, Flags )
    ! Set elements of Flags that correspond to nonzero elements of column C
    ! of Sparse to be true.  Flags is not initially made false -- see
    ! Sparse_Clear_Col_Flags.
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
    integer, intent(in) :: C               ! Which column
    logical, intent(inout) :: Flags(:)     ! The vector
    integer :: J
    j = sparse%cols(c)     ! Last element in column c
    if ( j == 0 ) return
    do
      j = sparse%e(j)%nc   ! Element in next column of column c
      flags(sparse%e(j)%r) = .true.
      if ( j == sparse%cols(c) ) exit
    end do
  end subroutine Sparse_Get_Col_Flags

  ! ------------------------------------------------  Vec_Dot_Col  -----
!  pure &
  pure real(rp) function Vec_Dot_Col ( Vector, Sparse, C ) result ( D )
    ! Compute dot product of Vector with column C of Sparse
    real(rp), intent(in) :: Vector(:)      ! The vector
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
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

  ! ---------------------------------------------  Vec_Dot_Sparse  -----
!  pure &
  subroutine Vec_Dot_Sparse ( Vector, Sparse, Prod )
    ! Multiply Sparse by Vector producing Product
    real(rp), intent(in) :: Vector(:)      ! The vector Sparse is to multiply
    class(sparse_t), intent(in) :: Sparse  ! The sparse matrix
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
! Revision 2.23  2020/05/07 21:26:13  vsnyder
! Output colon if requested in dump
!
! Revision 2.22  2018/12/13 01:15:27  vsnyder
! Add Rows argument to Row_Dot_Matrix
!
! Revision 2.21  2018/12/13 00:26:26  vsnyder
! Add Row_Dot_Matrix
!
! Revision 2.20  2018/11/16 22:52:46  vsnyder
! Don't add zero elements
!
! Revision 2.19  2018/11/07 21:24:09  vsnyder
! Checked it in last time before saving some changes to comments
!
! Revision 2.18  2018/11/07 21:22:32  vsnyder
! Add Allocate_Test and Deallocate_Test for Sparse_Element_t as private
! procedures.  Use them throughout.  Add "resize" option to Empty_Sparse.
!
! Revision 2.17  2018/10/27 01:35:57  vsnyder
! Add ExpLog and LogOnly arguments to dot products
!
! Revision 2.16  2018/10/25 19:00:51  vsnyder
! Spiff a dump
!
! Revision 2.15  2018/10/24 19:53:32  vsnyder
! Faster matrix-vector product
!
! Revision 2.14  2018/10/11 01:01:42  vsnyder
! More blank stuff to make transpose and non-transpose dumps more alike
!
! Revision 2.13  2018/10/11 00:43:29  vsnyder
! Inadvertently omitted a space in Dump_Sparse
!
! Revision 2.12  2018/10/11 00:33:31  vsnyder
! Add One_D_Col option to print 1D column subscript in Dump_Sparse
!
! Revision 2.11  2018/09/05 20:59:32  vsnyder
! Rename Row_Dot_Vec to Row_Dot_Vec_1D.  Add Row_Dot_Vec_2D and Row_Dot_Vec
! generic.  Rename Sparse_Dot_Vec to Sparse_Dot_Vec_1D.  Add Sparse_Dot_Vec_2D
! and Sparse_Dot_Vec generic.  Add Details argument ot List_Sparse.
!
! Revision 2.10  2018/08/21 01:52:14  vsnyder
! Add Exclude argument to sparse dumps
!
! Revision 2.9  2018/08/16 02:16:41  vsnyder
! Announce error on attempt to dump un-created Sparse_t
!
! Revision 2.8  2018/08/14 23:35:31  vsnyder
! Add Densify
!
! Revision 2.7  2018/08/03 23:17:49  vsnyder
! Make Sparse_Element_t public
!
! Revision 2.6  2018/05/24 03:20:26  vsnyder
! Respect the Width argument in the non-transpose dump case.  Do some
! cannonball polishing.
!
! Revision 2.5  2018/05/14 23:25:29  vsnyder
! Change to sparse eta representation
!
! Revision 2.4  2018/04/10 23:27:23  vsnyder
! Add Sparse_Clear_Col_And_Flags, Sparse_Get_Col_And_Flags.  Deallocate old
! Sparse%E before attempting to allocate one of a different size.
!
! Revision 2.3  2018/03/07 00:19:08  vsnyder
! Don't make names of procedures that are type-bound public.  Add
! Dump_Sparse_Array, Copy_Sparse, Invert_Column_Indices, Sparse_Clear_All_Flags,
! Sparse_Clear_Col, Sparse_Clear_Col_Flags, Sparse_Get_All_Flags, Sparse_Get_Col,
! Sparse_Get_Col_And_Sparsity, Sparse_Get_Col_Flags.  Add thumbnails.  Spiff
! the dumps.
!
! Revision 2.2  2017/11/28 23:58:02  vsnyder
! Add Empty_Sparse, List_Sparse, Row_Times_Vec, Sparse_Clear_Vec,
! Sparse_Dot_Matrix, Check_Column_Dup routines.  Add NRows element. Add
! transposed dump.  Correct some bugs.  Some cannonball polishing.
!
! Revision 2.1  2017/11/01 18:52:11  vsnyder
! Initial commit
!
