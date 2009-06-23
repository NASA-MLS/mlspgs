module BandedMatrix_m

  implicit NONE

  private

  integer, parameter :: RK = kind(0.0)

  type BandedMatrix
    integer, pointer :: FR(:) => NULL()   ! First and
    integer, pointer :: LR(:) => NULL()   !  last nonzero row per column
    integer, pointer :: FC(:) => NULL()   ! First and
    integer, pointer :: LC(:) => NULL()   !  last nonzero column per row
    real(rk), pointer :: V(:,:) => NULL() ! Values
  end type BandedMatrix

  public :: BandedMatrix, ClearBandedMatrix, CreateBandedMatrix
  public :: DestroyBandedMatrix

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains

  ! ------------------------------------------  ClearBandedMatrix  -----
  subroutine ClearBandedMatrix ( Matrix )
    ! Replace nonzeros with zeros.  Assume more nonzeros per
    ! column than row
    type(bandedMatrix), intent(inout) :: Matrix
    integer :: J
    do j = 1, size(matrix%fr)
      matrix%v(matrix%fr(j):matrix%lr(j),j) = 0.0_rk
    end do
    matrix%fr = size(matrix%fc)
    matrix%lr = 0
    matrix%fc = size(matrix%fr)
    matrix%lc = 0
  end subroutine ClearBandedMatrix

  ! -----------------------------------------  CreateBandedMatrix  -----
  subroutine CreateBandedMatrix ( Matrix, Rows, Columns )
    ! Create an empty banded matrix representing zeros.
    ! Destroy existing fields if new ones need to be different sizes.
    use Allocate_Deallocate, only: Deallocate_Test, Allocate_Test
    type(bandedMatrix), intent(inout) :: Matrix
    integer, intent(in) :: Rows, Columns

    if ( associated(matrix%v) ) then
      if ( size(matrix%v,1) /= rows .or. size(matrix%v,2) /= columns ) then
        call deallocate_test ( matrix%fr, 'matrix%fr', moduleName )
        call deallocate_test ( matrix%lr, 'matrix%lr', moduleName )
        call deallocate_test ( matrix%fc, 'matrix%fc', moduleName )
        call deallocate_test ( matrix%lc, 'matrix%lc', moduleName )
        call deallocate_test ( matrix%v , 'matrix%v',  moduleName )
      end if
    end if
    if ( .not. associated(matrix%v) ) then
      call allocate_test ( matrix%fr, columns, 'matrix%fr', moduleName )
      call allocate_test ( matrix%lr, columns, 'matrix%lr', moduleName )
      call allocate_test ( matrix%fc, rows,    'matrix%fc', moduleName )
      call allocate_test ( matrix%lc, rows,    'matrix%lc', moduleName )
      call allocate_test ( matrix%v, rows, columns, 'matrix%v', moduleName )
    end if
    matrix%fr = rows
    matrix%lr = 0
    matrix%fc = columns
    matrix%lc = 0
  end subroutine CreateBandedMatrix

  ! ----------------------------------------  DestroyBandedMatrix  -----
  subroutine DestroyBandedMatrix ( Matrix )
    ! Deallocate the fields of a banded matrix.
    use Allocate_Deallocate, only: Deallocate_Test
    type(bandedMatrix), intent(inout) :: Matrix
    call deallocate_test ( matrix%fr, 'matrix%fr', moduleName )
    call deallocate_test ( matrix%lr, 'matrix%lr', moduleName )
    call deallocate_test ( matrix%fc, 'matrix%fc', moduleName )
    call deallocate_test ( matrix%lc, 'matrix%lc', moduleName )
    call deallocate_test ( matrix%v , 'matrix%v',  moduleName )
  end subroutine DestroyBandedMatrix

  ! ------------------------------------------------  Not_Used_Here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module BandedMatrix_m

! $Log$
! Revision 2.2  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2007/06/20 22:25:08  vsnyder
! Initial commit
!
