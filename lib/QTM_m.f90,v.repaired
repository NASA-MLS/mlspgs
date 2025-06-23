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
  public :: Stack_t, ZOT_t, ZOT_V_t
  public :: H_t, H_Geoc, H_Geod ! From Geolocation_0 module
  ! Most procedures are type bound.  The specific and generic
  ! bindings are public.

  ! Procedures:
  public :: Expand_ZOT, Geo_To_ZOT, Get_Octant, High_Bit_Index, Init_Stack, &
    & QTM_Decode

  ! Kind of geolocation variables (from Geolocation_0 module):
  public :: RG

  integer, parameter :: Bits = ceiling ( digits(0) * log10(radix(0)+0.0d0) / log10(2.0) )
  integer, parameter, public :: QTM_Depth = (bits-4)/2
  ! QTM_Depth = 13 for 31-bit integers gives 1220 meter edge on smallest facet

  ! Kind of QTM facet ID
  integer, parameter, public :: QK = &
    & selected_int_kind(floor(log10(2.0_rg) * (2 * QTM_depth + 3)))

  ! Convert 1 + ZOT coordinate to an integer
  real(rg), parameter, public :: Scale_ZOT_x = 2.0_rg**(QTM_depth-1)
  real(rg), parameter, public :: Scale_ZOT_y = 2.0_rg**(2*QTM_depth)
  ! Convert an integer to 1 + a ZOT coordinate
  real(rg), parameter, public :: UnScale_ZOT_x = 2.0_rg**(1-QTM_depth)
  real(rg), parameter, public :: UnScale_ZOT_y = 2.0_rg**(-QTM_depth)

  ! Zenithial Ortho Triangular (ZOT) Projection coordinates
  type :: ZOT_t
    real(rg) :: X, Y
  contains
    procedure :: Condense_ZOT
    procedure :: Get_Octant
    procedure :: QTM_Encode_Level
    procedure :: QTM_Encode_Tol
    generic :: QTM_Encode => QTM_Encode_Level
    generic :: QTM_Encode => QTM_Encode_Tol
    procedure :: ZOT_To_Geo
    procedure :: ZOT_eq
    generic :: Operator ( == ) => ZOT_eq
    procedure :: ZOT_ne
    generic :: Operator ( /= ) => ZOT_ne
  end type ZOT_t

  type, extends(ZOT_t) :: ZOT_V_t ! ZOT plus a vertical coordinate
    real(rg) :: V
  end type ZOT_V_t

  type :: Stack_t
    integer :: top(8:15) = 0 ! Top frame for each octant
    type(zot_t), dimension(3,qtm_depth,8:15) :: Z
    integer, dimension(qtm_depth,8:15) :: XNode, YNode
    integer(qk) :: QID(qtm_depth,8:15) = -1 ! -1 => stack needs initialization
  contains
    procedure :: C_ZOT ! get ZOT coordinates of centroid
    procedure :: P_ZOT ! get ZOT coordinates of pole node
    procedure :: X_ZOT ! get ZOT coordinates of X node
    procedure :: Y_ZOT ! get ZOT coordinates of Y node
  end type Stack_t

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  pure elemental &
  integer(qk) function Condense_ZOT ( Z )
    !{ ZOT coordinates are in the range $[-1,+1]$.  At a refinement level of
    !  $l$, there are $2^l$ equally spaced values of each coordinate, with a
    !  spacing of $2^{l-1}$.  Since a QTM ID is bounded by $2^{2l+3}$, the
    !  integer $2^{q-1}(1+x) + 2^{2q}(1+y)$, where $q=${\tt QTM_depth}, is an
    !  unique identifier for a ZOT coordinate that fits in one integer of
    !  kind QK.  To avoid aliases, use $|x|$ if $|y|=1$ and use $|y|$ if
    !  $|x|=1$.

    class(zot_t), intent(in) :: Z
    real(rg) :: X, Y ! Disambiguated ZOT coordinates in 0 .. 2
    ! Disambiguate ZOT coordinates and put them in 0 .. 2
    x = merge(z%x,abs(z%x),abs(z%y)/=1.0_rg) + 1.0_rg
    y = merge(z%y,abs(z%y),abs(z%x)/=1.0_rg) + 1.0_rg
    ! Convert disambiguated ZOT coordinates to integers, and combine them
    condense_ZOT = int( x * scale_zot_x ) + int ( y * scale_zot_y )
  end function Condense_ZOT 

  pure elemental &
  type(zot_t) function Expand_ZOT ( IZ ) result ( Z )
    integer(qk), intent(in) :: IZ
    integer(qk) :: IX, IY
    integer, parameter :: Mask = 2**QTM_depth - 1
    ix = iand(iz,mask)
    iy = shiftr(iz,QTM_depth)
    z = zot_t ( ix * unscale_ZOT_x - 1.0_rg, iy * unscale_ZOT_y - 1.0_rg )
  end function Expand_ZOT

  pure elemental &
  type(zot_t) function Geo_To_ZOT ( Geo ) result ( Z )
    ! Given longitude and latitude (degrees), compute ZOT coordinates.
    class(h_t), intent(in) :: Geo ! Longitude and latitude; it doesn't matter
                                  ! whether latitude is geocentric or geodetic
    real(rg) :: Lon, Dxy, Temp
    integer :: LonInt
    ! Assume |lat| <= 90.
    ! Get fractional colatitude
    dxy = 1.0 - abs(geo%lat) / 90.0
    ! Put longitude in [0, 360).
    lon = mod(geo%lon%d,360.0_rg)
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
  end function Geo_To_ZOT

  pure elemental &
  integer function Get_Octant ( Z ) result ( Oct )
    ! Return the octant 8..15 occupied by any Z, or an error code -1..-15.
    ! This is 7 more than the octant number as defined by Dutton, so that
    ! the high-order nonzero bit of a QID is always 3 bits to the left the
    ! low-order bit of the QID.
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

  pure elemental &
  integer function High_Bit_Index ( Q ) result ( H )
    integer(qk), intent(in) :: Q
    do h = 1, bit_size(q)
      if ( shiftr(q,h) == 0 ) return
    end do
  end function High_Bit_Index

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
    integer :: Nx, Ny, Oct, Oct1, Octn
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
      ny = ( oct - 8 ) / 4 + 2    ! 2 for northern hemisphere, 3 for southern
      nx = 5 - ny                 ! 3 for northern hemisphere, 2 for southern
      stack%xnode(1,oct) = nx     ! x-node basis number
      stack%ynode(1,oct) = ny     ! y-node basis number
      ! Now that nodes are numbered, install coordinates for them
      stack%z(6-nx-ny,1,oct) = zPoleNode(oct)
      stack%z(nx,1,oct) = zXNode(oct)
      stack%z(ny,1,oct) = zYNode(oct)
      ! Fill the remaining stack frames with values indicating uninitialized
      stack%qid(2:,oct) = 9
      stack%xnode(2:,oct) = 0
      stack%ynode(2:,oct) = 0
      stack%z(:,2:,oct) = zot_t ( 0.0, 0.0 )
    end do
    stack%top = 1
  end subroutine Init_Stack

  ! Get ZOT coordinates of centroid
  pure elemental &
  type(zot_t) function C_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = zot_t ( ( sum(s%z(:,sp,oct)%x) ) / 3.0, ( sum(s%z(:,sp,oct)%y) ) / 3.0 )
  end function C_ZOT

  ! Get ZOT coordinates of pole node
  pure elemental &
  type(zot_t) function P_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    integer :: NP
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    np = 6 - ( s%xNode(sp,oct) + s%yNode(sp,oct) )
    z = s%z(np,sp,oct)
  end function P_ZOT

  ! Get ZOT coordinates of X node
  pure elemental &
  type(zot_t) function X_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = s%z(s%xNode(sp,oct),sp,oct)
  end function X_ZOT

  ! Get ZOT coordinates of Y node
  pure elemental &
  type(zot_t) function Y_ZOT ( S, Oct, QL ) result ( Z )
    class(stack_t), intent(in) :: S
    integer, intent(in) :: Oct ! Octant
    integer, intent(in), optional :: QL  ! Stack pointer, default top
    integer :: SP ! Stack pointer, default TOP of QL is not present
    sp = s%top(oct)
    if ( present(ql) ) sp = ql
    z = s%z(s%yNode(sp,oct),sp,oct)
  end function Y_ZOT

  pure elemental &
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
      pn = stack%z(np,ql,oct)
      if ( h == 0 ) exit
      dz = 0.5 * dz
      ql = ql + 1
      h = h - 2
      node = iand(shiftr(qid,h),3)
      if ( sameLoc .and. ql <= stack%top(oct) .and. stack%qid(ql,oct) == node ) then
        if ( present(lastMatch) ) lastMatch = ql
      else
        sameLoc = .false.
        call MakeStackConsistent ( Stack, ql, nP, pn, nX, nY, oct, node )
      end if
      stack%top(oct) = ql
    end do
  end subroutine QTM_Decode

! pure & ! Cannot be pure because Stack is intent(inout)
  integer(qk) function QTM_Encode_Level ( Z, Level, Stack, Oct, LastMatch )
    ! Get the QTM ID corresponding to Z.  Its digits are put on the stack.
    ! The octant is the four bits of the result value that begin with the
    ! highest-order nonzero bit.
    ! The remaining digits are taken from the stack, bottom up, and
    ! occupy two bits in the result.
    ! If the result is negative, its value is the error indicator from
    ! Get_Octant
    class(zot_t), intent(in) :: Z ! Assume z%x and z%y in -1 ... +1
    integer, intent(in) :: Level  ! Tol is 2**(-Level)
    type(stack_t), intent(inout) :: Stack
    integer, intent(out), optional :: Oct       ! Octant containing Z
    integer, intent(out), optional :: LastMatch ! Stack depth of last match
    qtm_encode_level = z%QTM_Encode ( 2.0_rg**(1-level), stack, oct, lastMatch )
  end function QTM_Encode_Level

! pure & ! Cannot be pure because Stack is intent(inout)
  integer(qk) function QTM_Encode_Tol ( Z, Tol, Stack, Oct, LastMatch )
    ! Get the QTM ID corresponding to Z.  Its digits are put on the stack.
    ! The octant is the four bits of the result value that begin with the
    ! highest-order nonzero bit.
    ! The remaining digits are taken from the stack, bottom up, and
    ! occupy two bits in the result.
    ! If the result is negative, its value is the error indicator from
    ! Get_Octant
    class(zot_t), intent(in) :: Z ! Assume z%x and z%y in -1 ... +1
    real(rg), intent(in) :: Tol   ! Assume tol is in 2**(-QTM_Depth) ... 1
    type(stack_t), intent(inout) :: Stack
    integer, intent(out), optional :: Oct       ! Octant containing Z
    integer, intent(out), optional :: LastMatch ! Stack depth of last match
    real(rg) :: Dx, Dy, Dz
    type(zot_t) :: Pn
    integer :: My_Oct, Node, Np, Nx, Ny, QL
    logical :: SameLoc
    my_oct = z%get_octant()
    if ( present(oct) ) oct = my_oct
    if ( my_oct < 0 ) then
      ! Z is defective.  z%x or z%y is outside the range -1 ... +1
      qtm_encode_tol = my_oct
      return
    end if
    if ( my_oct /= stack%qid(1,my_oct) ) call init_stack ( stack, my_oct )
    ql = 1           ! Initial QTM level
    sameLoc = .true. ! Initial flag indicating stack validity
    dz = 1.0         ! Initial ZOT x,y edge length for octant
    do while ( dz > tol .and. QL < size(stack%qid,1) )
      ! Get IDs (1,2,3) of X- and Y- nodes (nx,ny) from level QL of stack
      nx = stack%xnode(ql,my_oct)
      ny = stack%ynode(ql,my_oct)
      np = 6 - ( nx + ny ) ! polenode ID
      ! Retrieve ZOT x,y of polenode from level QL of stack
      pn = stack%z(np,ql,my_oct)
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
      if ( sameLoc .and. ql <= stack%top(my_oct) .and. stack%qid(ql,my_oct) == node ) then
        if ( present(lastMatch) ) lastMatch = ql
      else
        sameLoc = .false.
        call MakeStackConsistent ( Stack, QL, nP, pn, nX, nY, my_oct, Node )
      end if
      stack%top(my_oct) = ql
    end do
    QTM_Encode_tol = my_oct
    do ql = 2, stack%top(my_oct)
      QTM_Encode_tol = 4 * QTM_Encode_tol + stack%qid(ql,my_oct)
    end do
  end function QTM_Encode_Tol

  pure elemental &
  logical function ZOT_eq ( Z1, Z2 )
    class(zot_t), intent(in) :: Z1
    class(zot_t), intent(in) :: Z2
    zot_eq = z1%x == z2%x .and. z1%y == z2%y
  end function ZOT_eq

  pure elemental &
  logical function ZOT_ne ( Z1, Z2 )
    class(zot_t), intent(in) :: Z1
    class(zot_t), intent(in) :: Z2
    zot_ne = z1%x /= z2%x .or. z1%y /= z2%y
  end function ZOT_ne

  pure elemental &
  function ZOT_To_Geo ( Z ) result ( Geo )
    ! Given ZOT coordinate, compute longitude and latitude (degrees)
    class(zot_t), intent(in) :: Z
    type(h_t) :: Geo
    ! It doesn't matter whether the latitude is geocentric or geodetic;
    ! the same value is produced.
    real(rg) :: dx, dxy, dy
    integer :: Oct
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
        geo%lon%d = 90.0 * dy / dxy
      else
        geo%lon%d = 90.0 * dx / dxy
      end if
    else
      if ( mod(oct,2) /= 0 ) then
        geo%lon%d = 90.0 * dx / dxy
      else
        geo%lon%d = 90.0 * dy / dxy
      end if
    end if
    ! Offset longitude based upon octet number
    geo%lon%d = geo%lon%d + mod((mod(oct,4)-1) * 90, 360)
    ! Negate west longitudes
    if ( geo%lon%d > 180.0 ) geo%lon%d = geo%lon%d - 360.0
  end function ZOT_To_Geo

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
    xn = stack%z(nx,ql-1,oct)
    yn = stack%z(ny,ql-1,oct)
    ! Copy Pn%[xy], Xn%[xy], Yn%[xy] to next level, overriding
    ! those values for the two nonselected nodes
    if ( node /= np ) then ! Pole shifts to hypotenuse median
      stack%z(np,ql,oct) = zot_t ( (xn%x + yn%x) / 2, (xn%y + yn%y) / 2 )
    else
      stack%z(np,ql,oct) = pn
    end if
    if ( node /= ny ) then ! Y node shifts to pole-x median
      stack%z(ny,ql,oct) = zot_t ( (xn%x + pn%x) / 2, pn%y )
    else
      stack%z(ny,ql,oct) = yn
    end if
    if ( node /= nx ) then ! X node shifts to pole-y median
      stack%z(nx,ql,oct) = zot_t ( pn%x, (yn%y + pn%y) / 2 )
    else
      stack%z(nx,ql,oct) = xn
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
! Revision 2.6  2016/09/23 01:54:34  vsnyder
! Move dumps to QTM_Dumps
!
! Revision 2.5  2016/03/25 00:47:07  vsnyder
! Lon component now needs to acces its %d component
!
! Revision 2.4  2015/12/30 23:50:20  vsnyder
! Make Get_Octant public (not just type bound).  Add ZOT_V_t.  Change name
! of QTM_Encode to QTM_Encode_Tol and add QTM_Encode_Level.  Make procedures
! pure and elemental where possible, and add comments to explain why it's
! not possible in two cases.
!
! Revision 2.3  2015/12/07 20:13:07  vsnyder
! Simplify Dump_QID
!
! Revision 2.2  2015/12/01 02:56:00  vsnyder
! Add Dump_QID, Expand_ZOT, High_Bit_Index, Condense_ZOT, ZOT_eq, ZOT_ne.
! Change octant from 0..7 to 8..15 to make high bit index unambiguous.
! Change QTM_Depth from 14 to 13 because octant is four bits, not 3.
! Change stack's ZOT from ZP, ZX, ZY to Z(3) to correct errors in
! Make_Stack_Consistent.
!
! Revision 2.1  2015/11/13 19:45:12  vsnyder
! Initial commit
!
