! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module VectorsModule            ! Vectors in the MLS PGS suite
!=============================================================================

  ! This module provides the simple functionality for vector quantities in the
  ! MLS Level 2 software, and related programs.

  !---------------------------------------------------------------------------

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use STRING_TABLE, only: DISPLAY_STRING, GET_STRING, STRING_LENGTH

  implicit none
  public

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment (=)
    module procedure AssignVector
  end interface

  interface DUMP
    module procedure DUMP_VECTORS, DUMP_VECTOR_TEMPLATES
  end interface

  interface operator (+)
    module procedure AddVectors
  end interface

  interface operator (-)
    module procedure SubtractVectors
  end interface

  interface operator (*)
    module procedure ConstantXVector ! by a scalar on the left
  end interface

  interface operator ( .DOT. )
    module procedure DotVectors
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This type describes a vector template

  type VectorTemplate_T
     
    ! Administrative stuff
    integer :: Id              ! Id code for vector (for checking purposes)
    integer :: Name            ! Sub-rosa index of name, if any, else zero

    ! General information about the vector

    integer :: NoQuantities    ! Number of quantities in the vector
    integer :: TotalInstances  ! Number of horizontal instances in the vector
    integer :: TotalElements   ! Total of numbers of elements in the vector
    integer, dimension(:), pointer :: QUANTITIES => NULL() ! Indices of the
    !                            quantity templates in the quantities database
  end type VectorTemplate_T

  ! This type describes the subset of the values of a vector that
  ! correspond to a single quantity.

  type VectorValue_T
    type (QuantityTemplate_T), pointer :: TEMPLATE => NULL() ! Template for
    ! this quantity.
    real(r8), dimension(:,:), pointer :: VALUES => NULL() ! The dimensions of
    ! VALUES are Frequencies (or 1) * Vertical Coordinates (or 1), and
    ! Horizontal Instances (scan or profile or 1).  These are taken from
    ! (template%noChans * template%noSurfs, template%noInstances).
    integer, dimension(:,:), pointer :: MASK => NULL() ! MASK is used to
    ! control whether elements of vectors are of interest. If MASK is not
    ! associated, every element is of interest.  Otherwise,the dimensions of
    ! MASK are (size(values,1)+bit_size(mask)-1)/bit_size(mask) and
    ! size(values,2).  Bits of MASK are used to determine what is not
    ! interesting.  Zero means the corresponding element of VALUES is
    ! interesting, and one means it is not.
  end type VectorValue_T

  ! This type describes a vector.

  type Vector_T
    integer :: Name            ! Sub-rosa index of the vector name
    type (VectorTemplate_T), pointer :: TEMPLATE => NULL() ! In the template
    ! database
    type (VectorValue_T), dimension(:), pointer :: QUANTITIES => NULL() ! The
    ! dimension of QUANTITIES is the same as for the QUANTITIES field of the
    ! vector template.  Each element of QUANTITIES here corresponds to the
    ! one in the same position in the QUANTITIES field of the vector template.
  end type Vector_T

  ! This incrementing counter is used to set the id field for a vector template

  integer, save, private :: vectorTemplateCounter = 0

  integer, parameter, private :: B_sizer = 0
  integer, parameter, private :: B = bit_size(b_sizer) ! can't use bit_size(b)

  private :: CreateValues

contains ! =====     Public Procedures     =============================

  !-------------------------------------------------  AddToVector  -----
  subroutine AddToVector ( X, Y )  ! X = X + Y
    ! Dummy arguments:
    type(Vector_T), intent(inout) :: X
    type(Vector_T), intent(in) :: Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    do i = 1, size(x%quantities)
      x%quantities(i)%values = x%quantities(i)%values + y%quantities(i)%values
    end do
  end subroutine AddToVector

  !--------------------------------------------------  AddVectors  -----
  type (Vector_T) function AddVectors ( X, Y ) result (Z)
  ! Add two vectors, producing one having the same template (but no name,
  ! of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call CloneVector ( z, x )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = x%quantities(i)%values + y%quantities(i)%values
    end do
  end function AddVectors

  !---------------------------------  AddVectorTemplateToDatabase  -----
  integer function AddVectorTemplateToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector template to a database of such templates, 
  ! creating the database if necessary.

    ! Dummy arguments
    type (VectorTemplate_T), dimension(:), pointer :: DATABASE
    type (VectorTemplate_T), intent(in) :: ITEM

    ! Local variables
    type (VectorTemplate_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddVectorTemplateToDatabase = newSize
  end function AddVectorTemplateToDatabase

  !-----------------------------------------  AddVectorToDatabase  -----
  integer function AddVectorToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector template to a database of such templates, 
  ! creating the database if necessary.

    ! Dummy arguments
    type (Vector_T), dimension(:), pointer :: DATABASE
    type (Vector_T), intent(in) :: ITEM

    ! Local variables
    type (Vector_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddVectorToDatabase = newSize
  end function AddVectorToDatabase

  !------------------------------------------------  AssignVector  -----
  subroutine AssignVector ( Z, X )
  ! Destroy Z, then assign X to it, by using pointer assignment for the
  ! components.  DO NOT DO Z = Z!  Notice that CopyVector uses deep
  ! assignment. Notice that if we have a loop with Z = vector-expr, this
  ! destroys Z at each iteration, so it is necessary to call DestroyVectorInfo
  ! only after the loop, if at all.
    type(Vector_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: X
    call destroyVectorInfo ( z )
    z%name = x%name
    z%template => x%template
    z%quantities => x%quantities
  end subroutine AssignVector

  !--------------------------------------------------------  AXPY  -----
  type (Vector_T) function AXPY ( A, X, Y ) result (Z)
  ! Multiply the vector X by A and add Y to it, producing one having the
  ! same template (but no name, of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    real, intent(in) :: A
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call CloneVector ( z, x )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = &
        & a * x%quantities(i)%values + y%quantities(i)%values
    end do
  end function AXPY

  !---------------------------------------------------  ClearMask  -----
  subroutine ClearMask ( MASK, TO_CLEAR )
  ! Clear bits of MASK indexed by elements of TO_CLEAR.  If TO_CLEAR is
  ! absent, clear all of the bits of MASK.
    integer, intent(inout), dimension(:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_CLEAR
    integer :: I, W, P
    if ( present(to_clear) ) then
      do i = 1, size(to_clear)
        w = to_clear(i) / b
        p = mod(to_clear(i), b)
        mask(w+1) = ibclr(mask(w+1),p)
      end do
    else
      mask = 0
    end if
  end subroutine ClearMask

  !-------------------------------------------------  CloneVector  -----
  subroutine CloneVector ( Z, X )
  ! Destroy Z, except its name.
  ! Create the characteristics of a vector to be the same template as a
  ! given one (except it has no name).  Values are allocated, but not
  ! filled.  Z's mask is allocated if X's is allocated, but it is not filled.

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using Z after it is no
  ! longer needed. Otherwise, a memory leak will result.  Also see
  ! AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: X
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    call destroyVectorInfo ( z )
    z%template = x%template
    allocate ( z%quantities(size(x%quantities)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "z%quantities" )
    do i = 1, size(x%quantities)
      z%quantities(i)%template => x%quantities(i)%template
    end do
    call createValues ( z )
    do i = 1, size(x%quantities)
      if ( associated(x%quantities(i)%mask) ) &
        & call createMask ( z%quantities(i) )
    end do
  end subroutine  CloneVector

  !---------------------------------------------  ConstantXVector  -----
  type (Vector_T) function ConstantXVector ( A, X ) result (Z)
  ! Multiply the vector X by A, producing one having the same template
  ! (but no name, of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    real(r8), intent(in) :: A
    type(Vector_T), intent(in) :: X
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    call CloneVector ( z, x )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = a * x%quantities(i)%values
    end do
  end function ConstantXVector

  !-------------------------------------  ConstructVectorTemplate  -----
  subroutine ConstructVectorTemplate ( name, quantities, selected, &
    & vectorTemplate )

  ! This subroutine creates a vectorTemplate from a list of quantities.
  ! The default ordering is currently by quantity.  Later versions may
  ! have optional parameters to request other orderings.

    ! Dummy arguments
    integer, intent(in) :: NAME         ! Sub-rosa of vector template name
    type (QuantityTemplate_T), intent(in) :: quantities(:)
    integer, intent(in) :: selected(:)  ! Which quantities are selected?
    type (VectorTemplate_T), intent(out) :: vectorTemplate

    ! Local variables
    integer :: status

    ! Executable code
    vectorTemplate%name = name
    vectorTemplate%noQuantities = SIZE(selected)
    vectorTemplate%totalInstances = SUM(quantities(selected)%noInstances)
    vectorTemplate%totalElements = &
      & SUM(quantities(selected)%noInstances*quantities(selected)%instanceLen)
    
    ! Allocate some arrays

    allocate ( vectorTemplate%quantities(vectorTemplate%noQuantities), &
      & STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"Vector quantities" )

    ! Copy quantities over

    vectorTemplate%quantities = selected

    ! Increment the id counter and set the id field
    vectorTemplateCounter = vectorTemplateCounter + 1
    vectorTemplate%id = vectorTemplateCounter
  end subroutine ConstructVectorTemplate

  ! -------------------------------------------------  CopyVector  -----
  subroutine CopyVector ( Z, X, CLONE ) ! If CLONE is present and .true.,
  ! Destroy Z, deep Z = X, except the name of Z is not changed.  Otherwise,
  ! copy only the values and mask of X to Z

    type(Vector_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: X
    logical, intent(in), optional :: CLONE
    integer :: I
    logical MyClone
    myclone = .false.
    if ( present(clone) ) myclone = clone
    if ( myclone ) then
      call cloneVector ( Z, X )
    else
      if ( x%template%id /= z%template%id ) call MLSMessage &
        & ( MLSMSG_Error, ModuleName, 'Incompatible vectors in CopyVector' )
    end if
    do i = 1, size(x%quantities)
      z%quantities(i)%values = x%quantities(i)%values
      if ( associated (x%quantities(i)%mask ) ) &
        & z%quantities(i)%mask = x%quantities(i)%mask
    end do
  end subroutine CopyVector

  ! -------------------------------------------------  CreateMask  -----
  subroutine CreateMask ( VectorValue )
  ! Allocate the MASK array for a vector quantity.
    type(VectorValue_T), intent(inout) :: VectorValue
    call allocate_test ( vectorValue%mask, (size(vectorValue%values,1)+b-1)/b, &
      & size(vectorValue%values,2), "MASK in CreateMask", ModuleName )
    vectorValue%mask = 0 ! All vector elements are interesting
  end subroutine CreateMask

  ! -----------------------------------------------  CreateVector  -----
  type(Vector_T) function CreateVector &
    & ( vectorName, vectorTemplate, quantities ) result ( vector )

  ! This routine creates an empty vector according to a given template
  ! Its mask is not allocated.  Use CreateMask if one is needed.

    ! Dummy arguments
    integer, intent(in) :: vectorName   ! Sub_rosa index
    type (VectorTemplate_T), intent(in), target :: VectorTemplate ! For vector
    type (QuantityTemplate_T), dimension(:), intent(in), target :: Quantities

    ! Local variables
    integer :: QUANTITY                 ! Loop index
    integer :: STATUS                   ! From Allocate

    ! Executable code

    vector%name = vectorName
    vector%template => vectorTemplate
    allocate ( vector%quantities(vectorTemplate%noQuantities), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Vector quantities"  )
    do quantity = 1, vectorTemplate%noQuantities
      vector%quantities(quantity)%template => &
        & quantities(vectorTemplate%quantities(quantity))
    end do
    call createValues ( vector )
  end function CreateVector

  ! --------------------------------------  DestroyVectorDatabase  -----
  subroutine DestroyVectorDatabase ( database )

  ! This subroutine destroys a vector database

    ! Dummy argument
    type (Vector_T),  dimension(:), pointer :: database

    ! Local variables
    integer :: l2gpIndex, Status

    if ( associated(database) ) then
      do l2gpIndex = 1, SIZE(database)
        call DestroyVectorInfo(database(l2gpIndex))
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyVectorDatabase

  ! ------------------------------------------  DestroyVectorInfo  -----
  subroutine DestroyVectorInfo ( Vector )

  ! This routine destroys a vector created above

    ! Dummy arguments
    type (Vector_T), intent(inout) :: VECTOR

    ! Local Variables:
    integer :: STATUS

    ! Executable code

    vector%name = 0
    nullify ( vector%template )
    if ( .not. associated(vector%quantities) ) return
    call destroyVectorValue ( vector )
    deallocate ( vector%quantities, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_deallocate // "vector%quantities" )
  end subroutine DestroyVectorInfo

  ! ------------------------------------------  DestroyVectorMask  -----
  subroutine DestroyVectorMask ( Vector )

  ! This routine destroys the masks stored in the vector.

    ! Dummy arguments
    type (Vector_T), intent(inout) :: VECTOR

    ! Local Variables:
    integer :: I
    do i = 1, size(vector%quantities)
      call deallocate_test ( vector%quantities(i)%mask, &
        & "vector%quantities(i)%mask", ModuleName )
    end do
  end subroutine DestroyVectorMask

  ! ------------------------------  DestroyVectorTemplateDatabase  -----
  subroutine DestroyVectorTemplateDatabase ( database )

  ! This subroutine destroys a vector template database

    ! Dummy argument
    type (VectorTemplate_T), dimension(:), pointer :: database

    ! Local variables
    integer :: l2gpIndex, Status

    if ( associated(database) ) then
       do l2gpIndex = 1, SIZE(database)
          call DestroyVectorTemplateInfo ( database(l2gpIndex) )
       end do
       deallocate ( database, stat=status )
       if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyVectorTemplateDatabase

  !-----------------------------------  DestroyVectorTemplateInfo  -----
  subroutine DestroyVectorTemplateInfo ( VectorTemplate )

  ! This subroutine destroys a vector template created above

    ! Dummy arguments
    type (VectorTemplate_T), intent(inout) :: vectorTemplate

    ! Executable code

    call deallocate_test ( vectorTemplate%quantities, &
      & MLSMSG_deallocate // "vectorTemplate%quantities", ModuleName )

    vectorTemplate%noQuantities = 0
    vectorTemplate%totalInstances = 0
    vectorTemplate%totalElements = 0
    vectorTemplate%id = 0
    vectorTemplate%name = 0
  end subroutine DestroyVectorTemplateInfo

  !------------------------------------------  DestroyVectorValue  -----
  subroutine DestroyVectorValue ( Vector )
  ! Destroy the "values" field in all of the quantities in a vector.  This
  ! is useful when forming normal equations little-by-little.
    type(vector_T), intent(inout) :: Vector

    integer :: I

    do i = 1, size(vector%quantities)
      call deallocate_test ( vector%quantities(i)%values, &
        & "vector%quantities(i)%values", ModuleName )
    end do
  end subroutine DestroyVectorValue

  !--------------------------------------------------  DotVectors  -----
  real(r8) function DotVectors ( X, Y ) result (Z)
  ! Compute the inner product of two vectors.

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot .DOT. vectors having different templates" )
    z = 0.0_r8
    do i = 1, size(x%quantities)
      z = z + sum( x%quantities(i)%values * y%quantities(i)%values )
    end do
  end function DotVectors

  ! -----------------------------------------------  DUMP_VECTORS  -----
  subroutine DUMP_VECTORS ( VECTORS )
    type(Vector_T), intent(in) :: VECTORS(:)
    integer :: I, J
    call output ( 'VECTORS: SIZE = ' )
    call output ( size(vectors), advance='yes' )
    do i = 1, size(vectors)
      call output ( i, 4 )
      call output ( ': ' )
      if ( vectors(i)%name /= 0 ) then
        call output ( ' Name = ' )
        call display_string ( vectors(i)%name )
      end if
      if ( vectors(i)%template%name /= 0 ) then
        call output ( ' Template_Name = ' )
        call display_string ( vectors(i)%template%name )
      end if
      call output ( ' Template_ID = ' )
      call output ( vectors(i)%template%id, advance='yes' )
      do j = 1, size(vectors(i)%quantities)
        call output ( j, 4 )
        call output ( "~" )
        if ( vectors(i)%quantities(j)%template%name /= 0 ) then
          call output ( ' Quantity_Template_Name = ' )
          call display_string ( vectors(i)%quantities(j)%template%name )
        end if
        call output ( ' Quantity_Template_ID = ' )
        call output ( vectors(i)%quantities(j)%template%id, advance='yes' )
        call dump ( vectors(i)%quantities(j)%values, '      Elements = ' )
      end do ! j
    end do ! i
  end subroutine DUMP_VECTORS

  ! --------------------------------------  DUMP_VECTOR_TEMPLATES  -----
  subroutine DUMP_VECTOR_TEMPLATES ( VECTOR_TEMPLATES )
    type(VectorTemplate_T), intent(in) :: VECTOR_TEMPLATES(:)
    integer :: I
    call output ( 'VECTOR_TEMPLATES: SIZE = ' )
    call output ( size(vector_templates), advance='yes' )
    do i = 1, size(vector_templates)
      call output ( i, 4 )
      call output ( ': Id = ' )
      call output ( vector_templates(i)%id )
      if ( vector_templates(i)%name /= 0 ) then
        call output ( ' Name = ' )
        call display_string ( vector_templates(i)%name )
      end if
      call output ( ' NoQuantities = ' )
      call output ( vector_templates(i)%noQuantities )
      call output ( ' TotalInstances = ' )
      call output ( vector_templates(i)%totalInstances )
      call output ( ' TotalElements = ' )
      call output ( vector_templates(i)%totalElements, advance='yes' )
      call dump ( vector_templates(i)%quantities, '      Quantities = ' )
    end do
  end subroutine DUMP_VECTOR_TEMPLATES

  ! ------------------------------------------  GetVectorQuantity  -----
  function GetVectorQuantity ( vector, quantity )

  ! This function returns a pointer to the information about one quantity
  ! within a vector.

    ! Dummy arguments
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: Quantity                 ! Quantity index

    ! Result
    type(VectorValue_T), pointer :: GetVectorQuantity

    if ( quantity < 1 .or. quantity > size(vector%quantities) ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & '"Quantity" is out of range as a subscript for "vector%quantities".' )
    GetVectorQuantity => vector%quantities(quantity)

  end function GetVectorQuantity

  ! ------------------------------------  GetVectorQuantityByType  -----
  function GetVectorQuantityByType ( vector, quantityType, &
    & molecule, radiometer )

  ! Given a quantity type index (l_...), this function returns the first
  ! quantity within the vector that has that type.  If molecule and/or
  ! radiometer are supplied, the quantity that has the specified type, as
  ! well as the specified molecule and/or radiometer index, is returned.

    ! Dummy arguments
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: QuantityType ! Quantity type index (l_...)
    integer, intent(in), optional :: Molecule     ! Molecule index (l_...)
    integer, intent(in), optional :: Radiometer   ! Radiometer index
    ! Result
    type (VectorValue_T), pointer :: GetVectorQuantityByType

    GetVectorQuantityByType => vector%quantities( getVectorQuantityIndexByType &
      & ( vector, quantityType, molecule, radiometer ) )
  end function GetVectorQuantityByType

  ! -------------------------------  GetVectorQuantityIndexByName  -----
  integer function GetVectorQuantityIndexByName ( vector, quantityName )

  ! Given a quantity name's sub-rosa index, this function returns the index
  ! of the quantity within the vector that has that name.

    ! Dummy arguments
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: QuantityName ! Quantity name sub-rosa index

    ! Local variables
    character(len=127) :: MSG
    integer :: Search

    ! Executable code
    do search = 1, size(vector%quantities)
      if ( quantityName == vector%quantities(search)%template%name ) then
        GetVectorQuantityIndexByName = search
    return
      end if
    end do
    call get_string ( quantityName, msg )
    msg(string_length(quantityName)+2:) = 'is not a quantity in vector'
    call get_string ( vector%name, msg(len_trim(msg)+2:) )
    call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )

  end function GetVectorQuantityIndexByName

  ! -------------------------------  GetVectorQuantityIndexByType  -----
  integer function GetVectorQuantityIndexByType ( vector, quantityType, &
    & molecule, radiometer )

  ! Given a quantity type index (l_...), this function returns the index
  ! of the first quantity within the vector that has that type.  If
  ! molecule and/or radiometer are supplied, the quantity that has the
  ! specified type, as well as the specified molecule and/or radiometer
  ! index, is returned.

    ! Dummy arguments
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: QuantityType ! Quantity type index (l_...)
    integer, intent(in), optional :: Molecule     ! Molecule index (l_...)
    integer, intent(in), optional :: Radiometer   ! Radiometer index

    ! Local variables
    character(len=127) :: MSG
    integer :: Search

    ! Executable code
    do search = 1, size(vector%quantities)
      if ( quantityType == vector%quantities(search)%template%quantityType ) then
        if ( present(molecule) ) then
          if ( vector%quantities(search)%template%molecule /= molecule ) cycle
        end if
        if ( present(radiometer) ) then
          if ( vector%quantities(search)%template%radiometerIndex /= &
            &  radiometer ) cycle
        end if
        GetVectorQuantityIndexByType = search
    return
      end if
    end do
    msg = 'There is no quantity in vector '
    if ( vector%name /= 0 ) then
      call get_string ( vector%name, msg(len_trim(msg)+2:) )
    else
      msg(len_trim(msg)+2:) = '[unnamed]'
    end if
    msg = trim(msg) // ' that has the required type'
    call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )

  end function GetVectorQuantityIndexByType

  !---------------------------------------------  MultiplyVectors  -----
  subroutine MultiplyVectors ( X, Y, Z )
  ! If Z is present, destroy Z and clone a new one from X, then
  ! Z = X # Y where # means "element-by-element"; otherwise X = X # Y

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    type(Vector_T), intent(in) :: Y
    type(Vector_T), intent(out), optional, target :: Z
    ! Local Variables:
    integer :: I                        ! Subscript and loop inductor
    type(Vector_T), pointer :: Result   ! associated to either X or Z
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot multiply vectors having different templates" )
    if ( present(z) ) then
      call CloneVector ( z, x )
      result => z
    else
      result => x
    end if
    do i = 1, size(x%quantities)
      result%quantities(i)%values = &
        & x%quantities(i)%values * y%quantities(i)%values
    end do
  end subroutine MultiplyVectors

  !-------------------------------------------------  ScaleVector  -----
  subroutine ScaleVector ( X, A, Y )
  ! Y = A*X if Y is present, else X = A*X.

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    real(r8), intent(in) :: A
    type(Vector_T), intent(out), optional, target :: Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    type(Vector_T), pointer :: Z
    ! Executable statements:
    z => x
    if ( present(y) ) then
      if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Scaled vector has different template from original' )
      z => y
    end if
    do i = 1, size(x%quantities)
      z%quantities(i)%values = a * x%quantities(i)%values
    end do
  end subroutine ScaleVector

  !-----------------------------------------------------  SetMask  -----
  subroutine SetMask ( MASK, TO_SET )
  ! Set bits of MASK indexed by elements of TO_SET.  If TO_SET is absent,
  ! set all of the bits of MASK.
    integer, intent(inout), dimension(:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_SET
    integer :: I, W, P
    if ( present(to_set) ) then
      do i = 1, size(to_set)
        w = to_set(i) / b
        p = mod(to_set(i), b)
        mask(w+1) = ibset(mask(w+1),p)
      end do
    else
      mask = not(0)
    end if
  end subroutine SetMask

  !------------------------------------------  SubtractFromVector  -----
  subroutine SubtractFromVector ( X, Y, Quant, Inst ) ! X = X - Y.

    ! Dummy arguments:
    type(Vector_T), intent(inout) :: X
    type(Vector_T), intent(in) :: Y
    integer, intent(in), optional :: Quant, Inst  ! If Quant\ is present,
    !  only that quantity is subtracted.  If furthermore Inst\ is present,
    !  only that instance is subtracted.  If Inst\ is present but Quant\
    !  is not, the entire vector is subtracted.

    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot subtract vectors having different templates" )
    if ( present(quant) ) then
      if ( present(inst) ) then
        x%quantities(quant)%values(:,inst) = &
          & x%quantities(quant)%values(:,inst) - &
          & y%quantities(quant)%values(:,inst)
      else
        x%quantities(quant)%values = x%quantities(quant)%values - &
          &                          y%quantities(quant)%values
      end if
    else
      do i = 1, size(x%quantities)
        x%quantities(i)%values = x%quantities(i)%values - &
          &                      y%quantities(i)%values
      end do
    end if
  end subroutine SubtractFromVector

  !---------------------------------------------  SubtractVectors  -----
  type (Vector_T) function SubtractVectors ( X, Y ) result (Z)
  ! Subtract Y from X, producing one having the same template (but no name,
  ! of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.  Also see AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot subtract vectors having different templates" )
    call CloneVector ( z, x )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = x%quantities(i)%values - y%quantities(i)%values
    end do
  end function SubtractVectors

! =====     Private Procedures     =====================================
  subroutine CreateValues ( Vector )
  ! Allocate space for the values of a vector.
    type(Vector_T), intent(inout) :: Vector
    integer :: QTY
    do qty = 1, size(vector%quantities)
      call allocate_test ( vector%quantities(qty)%values, &
        & vector%quantities(qty)%template%noChans * &
        & vector%quantities(qty)%template%noSurfs, &
        & vector%quantities(qty)%template%noInstances, &
        & "Vector%quantities(qty)%values", ModuleName )
    end do
  end subroutine
!=======================================================================
end module VectorsModule
!=======================================================================

!
! $Log$
! Revision 2.9  2001/01/26 19:00:02  vsnyder
! Periodic commit
!
! Revision 2.8  2001/01/19 23:49:59  vsnyder
! Periodic commit
!
! Revision 2.7  2001/01/10 21:03:14  vsnyder
! Periodic commit
!
! Revision 2.6  2001/01/03 02:01:30  vsnyder
! Add molecule/radiometer functionality to GetVectorQuantityIndexByType
!
! Revision 2.5  2000/12/04 23:43:59  vsnyder
! Move more of addItemToDatabase into the include
!
! Revision 2.4  2000/11/23 01:10:27  vsnyder
! Add "mask" field to specify columns to ignore when vector is row- or
! column-specifier for a matrix.
!
! Revision 2.3  2000/11/15 01:33:58  vsnyder
! Added copyVector, assignment(=)
!
! Revision 2.2  2000/11/10 00:24:24  vsnyder
! Changed VectorValue_t%values from rank-3 to rank-2
!
! Revision 2.1  2000/10/13 00:00:37  vsnyder
! Moved from mlspgs/l2 to mlspgs/lib
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

