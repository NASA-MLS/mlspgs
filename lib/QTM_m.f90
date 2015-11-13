! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module QTM_m
!=============================================================================

  ! Types and procedures for Quaternary Triangular Meshes, based upon
  ! "A Hierarchical Coordinate System for Geoprocessing and Cartography"
  ! by Geoffrey H. Dutton, Springer Lecture Notes in Earth Sciences (1999)
  ! ISBN 3-540-64980-8.

  ! Some typos in the book are corrected.  In the parts of QTMdecode and
  ! QTMencode labeled "Is state stack still consistent at this level", the
  ! subscript for the stack should be QL-1, not zero.  Then in the three tests
  ! that compare "node" with "pn", "xn", or "yn", ".NE." should be ".EQ."

  use Geolocation_0, only: H_t, H_Geoc, H_Geod, RG

  implicit NONE
  private

  ! Types:
  public :: ZOT_t, Stack_t
  public :: H_t, H_Geoc, H_Geod ! From Geolocation_0 module
  ! Most procedures are type bound.  The specific and generic
  ! bindings are public.

  ! Procedures:
  public :: Dump, Dump_Stack, Geo_To_Zot, Init_Stack, QTM_Decode

  interface Dump
    module procedure Dump_Stack
  end interface

  ! Kind of geolocation variables (from Geolocation_0 module):
  public :: RG

  integer, parameter :: Bits = ceiling ( digits(0) * log10(radix(0)+0.0d0) / log10(2.0) )
  integer, parameter, public :: QTM_Depth = (bits-3)/2
  ! QTM_Depth = 14 for 31-bit integers gives 610 meter edge on smallest facet

  ! Kind of QTM facet ID
  integer, parameter, public :: QK = &
    & selected_int_kind(floor(log10(2.0_rg) * (2 * QTM_depth + 3)))

  ! Zenithial Ortho Triangular (ZOT) Projection coordinates
  type :: ZOT_t
    real(rg) :: X, Y
  contains
    procedure :: Get_Octant
    procedure :: QTM_Encode
    procedure :: Zot_To_Geo
  end type ZOT_t

  type :: Stack_t
    integer :: top(8:15) = 0 ! Top frame for each octant
    type(zot_t), dimension(qtm_depth,8:15) :: ZP, ZX, ZY
    integer, dimension(qtm_depth,8:15) :: XNode, YNode
    integer(qk) :: QID(qtm_depth,8:15) = -1 ! -1 => stack needs to be initialized
  contains
    procedure :: C_ZOT ! get ZOT coordinates of centroid
    procedure :: P_ZOT ! get ZOT coordinates of pole node
    procedure :: X_ZOT ! get ZOT coordinates of X node
    procedure :: Y_ZOT ! get ZOT coordinates of Y node
  end type Stack_t

  character(len=*), parameter, public :: StackDumpHead = &
    & ' QL  pNode X         Y      xNode X         Y      yNode X         Y       QID'

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Dump_Stack ( Stack, Oct, Top, Advance, Head )
    use Output_m, only: Output
    type(stack_t), intent(in) :: Stack
    integer, intent(in), optional :: Oct
    logical, intent(in), optional :: Top  ! Dump only the top frame
    character(*), intent(in), optional :: Advance
    character(*), intent(in), optional :: Head ! Print the heading followed
                                               ! by trim(head)
    integer :: I, I1, L, O, Oct1, Octn, Pn, Xn, Yn
    character(127) :: Line
    character(3) :: MyAdv
    logical :: MyTop
    oct1 = 8
    octn = 15
    if ( present(oct) ) then
      oct1 = oct
      octn = oct
    end if
    myTop = .false. ! Dump the whole stack
    if ( present(top) ) myTop = top
    do o = oct1, octn
      i1 = 1
      if ( myTop ) i1 = stack%top(o)
      if ( present(head) ) then
        call output ( stackDumpHead )
        call output ( trim(head), advance=advance )
      end if
      do i = i1, stack%top(o)
        xn = stack%xNode(i,o)
        yn = stack%yNode(i,o)
        pn = 6 - ( xn + yn )
        write ( line, '(i3,4(i3,2f10.5),i3)' ) i, &
          & pn, stack%zp(i,o), &
          & xn, stack%zx(i,o), &
          & yn, stack%zy(i,o)
        l = len_trim(line) + 4
        if ( myTop ) then
          write ( line(l:), '(99i0)' ) stack%qid(1,o)-7, &
            & stack%qid(2:stack%top(o),o)
        else
          write ( line(l:), '(i0)' ) stack%qid(i,o) - merge(7,0,i==1)
        end if
        myAdv = 'yes'
        if ( i == stack%top(o) ) myAdv = advance
        call output ( trim(line), myAdv )
      end do
    end do
  end subroutine Dump_Stack

  pure &
  type(zot_t) function Geo_To_Zot ( Geo ) result ( Z )
    ! Given longitude and latitude (degrees), compute ZOT coordinates.
    class(h_t), intent(in) :: Geo ! Longitude and latitude; doesn't matter
                                  ! whether latitude is geocentric or geodetic
    real(rg) :: Lon, Dxy, Temp
    integer :: LonInt
    ! Assume |lat| <= 90.
    ! Get fractional colatitude
    dxy = 1.0 - abs(geo%lat) / 90.0
    ! Put longitude in [0, 360).
    lon = mod(geo%lon,360.0_rg)
    if ( lon < 0.0 ) lon = lon + 360.0
    lonInt = lon
    z%x = mod(lonInt,90) + lon - lonInt
    ! Compute L1 "Manhattan" distance
    z%x = z%x * ( dxy / 90.0 )
    z%y = dxy - z%x
    if ( geo%lat < 0.0 ) then ! Reverse direction in southern hemisphere
      temp = 1.0 - z%x
      z%x = 1.0 - z%y
      z%y = temp
    end if
    if ( mod(lonInt/90,2) /= 0 ) then ! Swap X & Y in octants 2,4,6,8
      temp = z%x
      z%x = z%y
      z%y = temp
    end if
    ! Negate Y value in top half of ZOT space
    if ( lon > 90.0 .and. lon <= 270.0 ) z%y = -z%y
    ! Negate X value in left half of ZOT space
    if ( lon > 180.0 ) z%x = -z%x
  end function Geo_To_Zot

  pure &
  integer function Get_Octant ( Z ) result ( Oct )
    ! Return the octant 7..15 occupied by any Z, or an error code -1..-15
    class(zot_t), intent(in) :: Z
    integer :: Error
    error = 0
    oct = 0
    if ( z%x > 0.0 ) then
      if ( z%x > 1.0 ) error = error - 1
      if ( z%y > 0.0 ) then
        if ( z%y > 1.0 ) error = error - 8
        oct = 8
      else
        if ( z%y < -1.0 ) error = error - 4
        oct = 9
      end if
    else
      if ( z%x < -1.0 ) error = error - 2
      if ( z%y < 0.0 ) then
        if ( z%y < -1.0 ) error = error - 4
        oct = 10
      else
        if ( z%y > 1.0 ) error = error - 8
        oct = 11
      end if
    end if
    if ( error /= 0 ) then
      oct = error
    else if ( abs(z%x) + abs(z%y) > 1.0 ) then
      ! In southern hemisphere, add 4 to octant number
      oct = oct + 4
    end if
  end function Get_Octant

  pure &
  subroutine Init_Stack ( Stack, Octant )
    type(stack_t), intent(out) :: Stack
    integer, intent(in), optional :: Octant
    type(zot_t), parameter :: zPoleNode(8:15) = &
      & [ zot_t( 0.0, 0.0), zot_t( 0.0, 0.0), zot_t( 0.0, 0.0), zot_t( 0.0, 0.0), &
      &   zot_t( 1.0, 1.0), zot_t( 1.0,-1.0), zot_t(-1.0,-1.0), zot_t(-1.0, 1.0) ]
    type(zot_t), parameter :: zXNode(8:15) = &
      & [ zot_t( 1.0, 0.0), zot_t( 1.0, 0.0), zot_t(-1.0, 0.0), zot_t(-1.0, 0.0), &
      &   zot_t( 0.0, 1.0), zot_t( 0.0,-1.0), zot_t( 0.0,-1.0), zot_t( 0.0, 1.0) ]
    type(zot_t), parameter :: zYNode(8:15) = &
      & [ zot_t( 0.0, 1.0), zot_t( 0.0,-1.0), zot_t( 0.0,-1.0), zot_t( 0.0, 1.0), &
      &   zot_t( 1.0, 0.0), zot_t( 1.0, 0.0), zot_t(-1.0, 0.0), zot_t(-1.0, 0.0) ]
    integer :: Ix, Iy, Oct, Oct1, Octn
    oct1 = 8
    octn = 15
    if ( present(octant) ) then
      oct1 = octant
      octn = octant
    end if
    do oct = oct1, octn
      ! Fill the first stack frame of the octant
      stack%qid(1,oct) = oct ! First QTM digit is octant number
      ! Octant polenodes' IDs are always == 1, so do not need to be stored
      ! We do need, however, to store the pole coordinates
      iy = ( oct - 8 ) / 4 + 1    ! 1 for northern hemisphere, 2 for southern
      ix = 3 - iy                 ! 2 for northern hemisphere, 1 for southern
      stack%xnode(1,oct) = ix + 1 ! x-node basis number
      stack%ynode(1,oct) = iy + 1 ! y-node basis number
      ! Now that nodes are numbered, install coordinates for them
      stack%zp(1,oct)  = zPoleNode(oct)
      stack%zx(1,oct) = zXnode(oct)
      stack%zy(1,oct) = zYnode(oct)
      ! Fill the remaining stack frames with values indicating uninitialized
      stack%qid(2:,oct) = 9
      stack%xnode(2:,oct) = 0
      stack%ynode(2:,oct) = 0
      stack%zp(2:,oct) = zot_t ( 0.0, 0.0 )
      stack%zx(2:,oct) = zot_t ( 0.0, 0.0 )
      stack%zy(2:,oct) = zot_t ( 0.0, 0.0 )
    end do
    stack%top = 1
  end subroutine Init_Stack

  ! Get ZOT coordinates of centroid
  pure &
  type(zot_t) function C_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = zot_t ( ( s%zp(sp,oct)%x + s%zx(sp,oct)%x + s%zy(sp,oct)%x ) / 3.0, &
              & ( s%zp(sp,oct)%y + s%zx(sp,oct)%y + s%zy(sp,oct)%y ) / 3.0 )
  end function C_ZOT

  ! Get ZOT coordinates of pole node
  pure &
  type(zot_t) function P_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = s%zp(sp,oct)
  end function P_ZOT

  ! Get ZOT coordinates of X node
  pure &
  type(zot_t) function X_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = s%zx(sp,oct)
  end function X_ZOT

  ! Get ZOT coordinates of Y node
  pure &
  type(zot_t) function Y_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = s%zy(sp,oct)
  end function Y_ZOT

  pure &
  subroutine QTM_Decode ( QID, Stack, LastMatch )
    ! Compute ZOT coordinates of a specified node of the facet identified by QID
    ! and leave them in the top frame of the stack.
    ! The highest-order one bit of the QID indicates the first bit of the octant
    ! number, because octants are numbered 8..15.  The difference between the
    ! first bit of the octant number and the number of bits determines the
    ! tolerance.
    integer(qk), intent(in) :: QID ! QTM ID to decode
    type(stack_t), intent(inout) :: Stack
    integer, intent(out), optional :: LastMatch ! Stack depth of last match
    real(rg) :: Dz
    integer :: H          ! High-order bit index, to find octant
    integer :: Node       ! Digit of QID
    integer :: Oct        ! Octant, from QID
    integer :: QL         ! Stack pointer
    integer :: Np, Nx, Ny ! Node IDs
    type(zot_t) :: Pn     ! Pole node
    logical :: SameLoc

    h = high_bit_index(QID) - 4
    oct = shiftr(QID,h)
    if ( oct /= stack%qid(1,oct) ) call init_stack ( stack, oct )
    ql = 1           ! Initial QTM level
    sameLoc = .true. ! Initial flag indicating stack validity
    dz = 1.0
    do
      ! Get IDs (1,2,3) of X- and Y- nodes (nx,ny) from level QL of stack
      nx = stack%xnode(ql,oct)
      ny = stack%ynode(ql,oct)
      np = 6 - ( nx + ny ) ! polenode ID
      ! Retrieve ZOT x,y of polenode from level QL of stack
      pn = stack%zp(ql,oct)
      if ( h == 0 ) exit
      dz = 0.5 * dz
      ql = ql + 1
      h = h - 2
      node = iand(shiftr(qid,h),3)
      if ( sameLoc .and. ql <= stack%top(oct) .and. stack%qid(ql,oct) == node ) then
        if ( present(lastMatch) ) lastMatch = ql
      else
        sameLoc = .false.
        call MakeStackConsistent ( Stack, QL, nP, pn, nX, nY, Oct, Node )
      end if
      stack%top(oct) = ql
    end do
  contains
    pure integer function High_Bit_Index ( Q ) result ( H )
      integer(qk), intent(in) :: Q
      do h = 1, bit_size(q)
        if ( shiftr(q,h) == 0 ) return
      end do
    end function High_Bit_Index
  end subroutine QTM_Decode

  integer(qk) function QTM_Encode ( Z, Tol, Oct, Stack, LastMatch )
    ! Get the QTM ID corresponding to Z.  Its digits are put on the stack.
    ! The octant is the four bits of the result value that begin with the
    ! highest-order nonzero bit.
    ! The remaining digits are taken from the stack, bottom up, and
    ! occupy two bits in the result.
    ! If the result is negative, its value is the error indicator from
    ! Get_Octant
    class(zot_t), intent(in) :: Z ! Assume z%x and z%y in -1 ... +1
    real(rg), intent(in) :: Tol   ! Assume tol is in 2**(-QTM_Depth) ... 1
    integer, intent(out) :: Oct   ! Octant containing Z
    type(stack_t), intent(inout) :: Stack
    integer, intent(out), optional :: LastMatch ! Stack depth of last match
    real(rg) :: Dx, Dy, Dz
    type(zot_t) :: Pn
    integer :: Node, Np, Nx, Ny, QL
    logical :: SameLoc
    oct = z%get_octant()
    if ( oct < 0 ) then
      ! Z is defective.  z%x or z%y is outside the range -1 ... +1
      qtm_encode = oct
      return
    end if
    if ( oct /= stack%qid(1,oct) ) call init_stack ( stack, oct )
    ql = 1           ! Initial QTM level
    sameLoc = .true. ! Initial flag indicating stack validity
    dz = 1.0         ! Initial ZOT x,y edge length for octant
    do while ( dz > tol .and. QL < size(stack%qid,1) )
      ! Get IDs (1,2,3) of X- and Y- nodes (nx,ny) from level QL of stack
      nx = stack%xnode(ql,oct)
      ny = stack%ynode(ql,oct)
      np = 6 - ( nx + ny ) ! polenode ID
      ! Retrieve ZOT x,y of polenode from level QL of stack
      pn = stack%zp(ql,oct)
      dz = 0.5 * dz ! Halve closeness criterion
      ! Compute displacement of Z from polenode
      dx = abs(z%x - pn%x)
      dy = abs(z%y - pn%y)
      ! Identify closest node to Z. Node 0 represents central facet.
      node = 0
      if ( dx + dy <= dz ) then
        node = np
      else if ( dx >= dz ) then
        node = nx
      else if ( dy >= dz ) then
        node = ny
      end if
      ql = ql + 1
      ! Is stack state consistent at this level?
      if ( sameLoc .and. ql <= stack%top(oct) .and. stack%qid(ql,oct) == node ) then
        if ( present(lastMatch) ) lastMatch = ql
      else
        sameLoc = .false.
        call MakeStackConsistent ( Stack, QL, nP, pn, nX, nY, Oct, Node )
      end if
      stack%top(oct) = ql
    end do
    QTM_Encode = oct
    do ql = 2, stack%top(oct)
      QTM_Encode = 4 * QTM_Encode + stack%qid(ql,oct)
    end do
  end function QTM_Encode

  pure &
  function Zot_To_Geo ( Z, Geodetic ) result ( Geo )
    ! Given ZOT coordinate, compute longitude and latitude (degrees)
    class(zot_t), intent(in) :: Z
    class(h_t), allocatable :: Geo
    logical, intent(in), optional :: Geodetic ! Default geocentric
    ! It doesn't matter whether the latitude is geocentric or geodetic;
    ! the same value is produced.  The difference is which type gets
    ! allocated for the result, i.e., what kind of latitude you want it
    ! to be called.  This is needed because function result types do not
    ! participate in generic resolution.
    real(rg) :: dx, dxy, dy
    integer :: Oct
    logical :: MyGeod
    myGeod = .false.
    if ( present(geodetic) ) myGeod = geodetic
    if ( myGeod ) then
      allocate ( h_geod :: geo )
    else
      allocate ( h_geoc :: geo )
    end if
    ! Assume |z%x| <= 1.0 and |z%y| <= 1.0
    oct = z%get_octant() - 7 ! 8..15 => 1..8
    dx = abs(z%x)
    dy = abs(z%y)
    if ( oct > 4 ) then ! Southern hemisphere
      dx = 1.0 - dx
      dy = 1.0 - dy
    end if
    dxy = dx + dy
    if ( dxy == 0.0 ) dxy = tiny(dxy) * 360.0 ! Avoid overflow near poles
    geo%lat = 90.0 * ( 1.0 - dxy )
    if ( oct > 4 ) then ! Southern hemisphere
      geo%lat = -geo%lat
      if ( mod(oct,2) /= 0 ) then
        geo%lon = 90.0 * dy / dxy
      else
        geo%lon = 90.0 * dx / dxy
      end if
    else
      if ( mod(oct,2) /= 0 ) then
        geo%lon = 90.0 * dx / dxy
      else
        geo%lon = 90.0 * dy / dxy
      end if
    end if
    ! Offset longitude based upon octet number
    geo%lon = geo%lon + mod((mod(oct,4)-1) * 90, 360)
    ! Negate west longitudes
    if ( geo%lon > 180.0 ) geo%lon = geo%lon - 360.0
  end function Zot_To_Geo

  pure &
  subroutine MakeStackConsistent ( Stack, QL, Np, Pn, Nx, Ny, Oct, Node )
    type(stack_t), intent(inout) :: Stack
    integer, intent(in) :: QL          ! Stack frame index
    integer, intent(in) :: Np          ! Current polenode index
    type(zot_t), intent(in) :: Pn      ! Polenode
    integer, intent(inout) :: Nx, Ny   ! Node IDs
    integer, intent(in) :: Oct         ! In which octant are we working?
    integer, intent(in) :: Node        ! QID to put into the stack if inconsistent
    integer :: Temp
    type(zot_t) :: Xn, Yn
    xn = stack%zx(ql-1,oct)
    yn = stack%zy(ql-1,oct)
    ! Copy Pn%[xy], Xn%[xy], Yn%[xy] to next level, overriding
    ! those values for the two nonselected nodes
    if ( node /= np ) then ! Pole shifts to hypotenuse median
      stack%zp(ql,oct) = zot_t ( (xn%x + yn%x) / 2, (xn%y + yn%y) / 2 )
    else
      stack%zp(ql,oct) = pn
    end if
    if ( node /= ny ) then ! Y node shifts to pole-x median
      stack%zy(ql,oct) = zot_t ( (xn%x + pn%x) / 2, pn%y )
    else
      stack%zy(ql,oct) = yn
    end if
    if ( node /= nx ) then ! X node shifts to pole-y median
      stack%zx(ql,oct) = zot_t ( pn%x, (yn%y + pn%y) / 2 )
    else
      stack%zx(ql,oct) = xn
    end if
    ! Renumber nodes for new level
    if ( node == nx ) then ! swap ny and np, but np is not needed again
      ny = np              ! Pole shifts to X node
    else if ( node == ny ) then ! swap nx and np, but np is not needed again
      nx = np              ! Pole shifts to Y node
    else if ( node == np ) then ! swap nx and ny
      temp = ny            ! Pole shifts to hypotenuse median
      ny = nx
      nx = temp
    end if
    stack%xnode(ql,oct) = nx
    stack%ynode(ql,oct) = ny
    stack%qid(ql,oct) = node
    stack%top(oct) = ql
  end subroutine MakeStackConsistent

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module QTM_m

! $Log$
! Revision 2.1  2015/11/13 19:45:12  vsnyder
! Initial commit
!
