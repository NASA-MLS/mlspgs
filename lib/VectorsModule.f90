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
  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER

  implicit none
  public

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment (=)
    module procedure AssignVector
  end interface

  interface DUMP
    module procedure DUMP_VECTOR, DUMP_VECTORS, DUMP_VECTOR_TEMPLATES
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
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This type describes a vector template

  type VectorTemplate_T
     
    ! Administrative stuff
    integer :: Id = 0          ! Id code for vector (for checking purposes)
    integer :: Name = 0        ! Sub-rosa index of name, if any, else zero

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
    type (QuantityTemplate_T) :: TEMPLATE ! Template for this quantity.
    integer :: index                    ! Index of this quantity into vector
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
    integer :: Name = 0        ! Sub-rosa index of the vector name
    type (VectorTemplate_T) :: TEMPLATE ! Template for this vector
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
    call CloneVector ( z, x, vectorNameText='_z' )
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
    z%template = x%template
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
    call CloneVector ( z, x, vectorNameText='_z' )
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
  subroutine CloneVector ( Z, X, VectorNameText )
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
    character(len=*), intent(in), optional :: VectorNameText
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    call destroyVectorInfo ( z )
    if ( present(vectorNameText) ) &
      & z%name = enter_terminal ( vectorNameText, t_identifier )
    z%template = x%template
    allocate ( z%quantities(size(x%quantities)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "z%quantities" )
    z%quantities%index = x%quantities%index
    do i = 1, size(x%quantities)
      z%quantities(i)%template = x%quantities(i)%template
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
    call CloneVector ( z, x, vectorNameText='_z' )
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
      call cloneVector ( Z, X, vectorNameText='_Z' )
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
    & ( vectorName, vectorTemplate, quantities, VectorNameText ) &
    & result ( vector )

  ! This routine creates an empty vector according to a given template
  ! Its mask is not allocated.  Use CreateMask if one is needed.

    ! Dummy arguments
    integer, intent(in) :: vectorName   ! Sub_rosa index
    type (VectorTemplate_T), intent(in), target :: VectorTemplate ! For vector
    type (QuantityTemplate_T), dimension(:), intent(in), target :: Quantities
    character(len=*), intent(in), optional :: VectorNameText

    ! Local variables
    integer :: QUANTITY                 ! Loop index
    integer :: STATUS                   ! From Allocate

    ! Executable code

    vector%name = vectorName
    if ( present(vectorNameText) ) &
      & vector%name = enter_terminal ( vectorNameText, t_identifier )
    vector%template = vectorTemplate
    allocate ( vector%quantities(vectorTemplate%noQuantities), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "Vector quantities"  )
    do quantity = 1, vectorTemplate%noQuantities
      vector%quantities(quantity)%index = quantity
      vector%quantities(quantity)%template = &
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
    nullify ( vector%template%quantities )
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

    if ( .not. associated(vector%quantities) ) return
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

  ! ------------------------------------------------  DUMP_VECTOR  -----
  subroutine DUMP_VECTOR ( VECTOR, DETAILS, NAME )
    type(Vector_T), intent(in) :: VECTOR
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump quantity values
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    integer :: J    ! Loop inductor, subscript
    integer :: MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( present(name) ) then
      call output ( name ); call output ( ', ' )
    end if
    if ( vector%name /= 0 ) then
      call output ( 'Name = ' )
      call display_string ( vector%name )
    end if
    if ( vector%template%name /= 0 ) then
      call output ( ' Template_Name = ' )
      call display_string ( vector%template%name )
    end if
    call output ( ' Template_ID = ' )
    call output ( vector%template%id, advance='yes' )
    do j = 1, size(vector%quantities)
      call output ( j, 4 )
      call output ( "~" )
      if ( vector%quantities(j)%template%name /= 0 ) then
        call output ( ' Qty_Template_Name = ' )
        call display_string ( vector%quantities(j)%template%name )
      end if
      call output ( ' Qty_Template_ID = ' )
      call output ( vector%quantities(j)%template%id )
      if ( myDetails > 0 ) then
        call dump ( vector%quantities(j)%values, ', Elements = ' )
        if ( associated(vector%quantities(j)%mask) ) then
          call dump ( vector%quantities(j)%mask, format='(z8)' )
        else
          call output ( '      Without mask', advance='yes' )
        end if
      else
        call output ( ', with' )
        if ( .not. associated(vector%quantities(j)%values) ) &
          & call output ( 'out' )
        call output ( ' values, with' )
        if ( .not. associated(vector%quantities(j)%mask ) ) &
          & call output ( 'out' )
        call output ( ' mask', advance='yes' )
      end if
      if ( associated(vector%quantities(j)%mask) ) then
        if ( myDetails > 0 ) then
          call dump ( vector%quantities(j)%mask, format='(z8)' )
        else
          call output ( '      With mask', advance='yes' )
        end if
      else
        call output ( '      Without mask', advance='yes' )
      end if
    end do ! j
  end subroutine DUMP_VECTOR

  ! -----------------------------------------------  DUMP_VECTORS  -----
  subroutine DUMP_VECTORS ( VECTORS, DETAILS, NAME )
    type(Vector_T), intent(in) :: VECTORS(:)
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump quantity values
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    integer :: I
    if ( size(vectors) > 1 ) then
      call output ( 'VECTORS: SIZE = ' )
      call output ( size(vectors), advance='yes' )
    end if
    do i = 1, size(vectors)
      call output ( i, 4 )
      call output ( ': ' )
      call dump_vector ( vectors(i), details, name )
    end do ! i
  end subroutine DUMP_VECTORS

  ! --------------------------------------  DUMP_VECTOR_TEMPLATES  -----
  subroutine DUMP_VECTOR_TEMPLATES ( VECTOR_TEMPLATES, DETAILS )
    type(VectorTemplate_T), intent(in) :: VECTOR_TEMPLATES(:)
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
    !                                        ! > 0  => Do dump arrays
    !                                        ! Default 1
    integer :: I, MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
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
      if ( myDetails > 0 ) &
        & call dump ( vector_templates(i)%quantities, '      Quantities = ' )
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
  function GetVectorQuantityByType ( vector, otherVector, quantityType, &
    & molecule, instrumentModule, radiometer, signal, &
    & sideband, foundInFirst, noError )

    ! Given a quantity type index (l_...), this function returns the first
    ! quantity within the vector that has that type.  If molecule and/or
    ! radiometer are supplied, the quantity that has the specified type, as
    ! well as the specified molecule and/or radiometer index, is returned.
    
    ! If otherVector is present it will look for the quantity there.

    ! Dummy arguments
    type (Vector_T), intent(in) :: VECTOR ! First vector to look in
    type (Vector_T), intent(in), optional :: OTHERVECTOR ! Second vector to look in
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in),  optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in),  optional :: INSTRUMENTMODULE ! Instrument module index
    integer, intent(in),  optional :: RADIOMETER   ! Radiometer index
    integer, intent(in),  optional :: SIGNAL       ! Signal index
    integer, intent(in),  optional :: SIDEBAND ! -1, 0, +1
    logical, intent(out), optional :: FOUNDINFIRST ! Set if found in first vector
    logical, intent(in),  optional :: NOERROR ! Don't give error if not found
    ! Result
    type (VectorValue_T), pointer :: GetVectorQuantityByType

    ! Local variable
    integer :: index
    logical :: myNoError

    myNoError = .false.
    if (present(noError)) myNoError = noError

    GetVectorQuantityByType => NULL()

    ! Look in the first vector
    index = GetVectorQuantityIndexByType ( vector, &
      & quantityType, molecule, instrumentModule, radiometer, signal, sideband, &
      &   noError = present(otherVector) .or. myNoError)
    if ( index /= 0 ) then
      if ( present (foundInFirst) ) foundInFirst = .true.
      GetVectorQuantityByType => vector%quantities(index)
    else
      ! Can only get here if not found in first vector and noError or other vector
      if ( present (otherVector) ) then
        index = GetVectorQuantityIndexByType ( otherVector, &
          &  quantityType, molecule, instrumentModule, radiometer, signal, sideband, &
          &  noError=myNoError )
        if ( present (foundInFirst) ) foundInFirst = .false.
        if ( index /= 0 ) GetVectorQuantityByType => otherVector%quantities( index )
      end if
    end if
  end function GetVectorQuantityByType

  ! ------------------------------- GetVectorQtyByTemplateIndex --i
  function GetVectorQtyByTemplateIndex(vector, quantityIndex, indexInVector )
    ! Given a vector and an index into the quantity templates, find quantity
    ! with matching template within vector.

    ! Dummy arguments
    type (vector_T), intent(in) :: vector
    integer, intent(in) :: quantityIndex
    integer, intent(out), optional :: indexInVector
    ! Result
    type (VectorValue_T), pointer :: GetVectorQtyByTemplateIndex

    ! Local variables
    integer :: i, myIndexInVector

    ! Executable code
    myIndexInVector=0
    GetVectorQtyByTemplateIndex => NULL()
    do i=1,vector%template%noQuantities
      if ( vector%template%quantities(i) == quantityIndex) then
        myIndexInVector=i
      end if
    end do
    if ( myIndexInVector /= 0 ) &
      & GetVectorQtyByTemplateIndex => &
      &   vector%quantities(myIndexInVector)
    if ( present ( indexInVector ) ) indexInVector = myIndexInVector
  end function GetVectorQtyByTemplateIndex

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
    & molecule, instrumentModule, radiometer, signal, sideband, noError )

  ! Given a quantity type index (l_...), this function returns the index
  ! of the first quantity within the vector that has that type.  If
  ! molecule and/or radiometer are supplied, the quantity that has the
  ! specified type, as well as the specified molecule and/or radiometer
  ! index, is returned.

    ! Dummy arguments
    type (Vector_T), intent(in) :: VECTOR
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in), optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in), optional :: INSTRUMENTMODULE ! Module index
    integer, intent(in), optional :: RADIOMETER   ! Radiometer index
    integer, intent(in), optional :: SIGNAL       ! Signal Index
    integer, intent(in), optional :: SIDEBAND ! -1, 0, +1
    logical, intent(in), optional :: NOERROR ! Don't give error if not found

    ! Local variables
    character(len=127) :: MSG
    integer :: SEARCH
    logical :: MYNOERROR

    myNoError = .false.
    if (present(noError)) myNoError = noError

    ! Executable code
    do search = 1, size(vector%quantities)
      if ( quantityType == vector%quantities(search)%template%quantityType ) then
        if ( present(molecule) ) then
          if ( vector%quantities(search)%template%molecule /= molecule ) cycle
        end if
        if ( present(instrumentModule) ) then
          if ( vector%quantities(search)%template%instrumentModule /= &
            &  instrumentModule ) cycle
        end if
        if ( present(radiometer) ) then
          if ( vector%quantities(search)%template%radiometer /= &
            &  radiometer ) cycle
        end if
        if ( present(signal) ) then
          if ( vector%quantities(search)%template%signal /= &
            &  signal ) cycle
        end if
        if ( present(sideband) ) then
          if ( vector%quantities(search)%template%sideband /= &
            &  sideband ) cycle
        end if
        GetVectorQuantityIndexByType = search
    return
      end if
    end do

    ! Not found, perhaps generate an error
    if (myNoError) then
      GetVectorQuantityIndexByType = 0
    else
      msg = 'There is no quantity in vector '
      if ( vector%name /= 0 ) then
        call get_string ( vector%name, msg(len_trim(msg)+2:) )
      else
        msg(len_trim(msg)+2:) = '[unnamed]'
      end if
      msg = trim(msg) // ' that has the required type'
      call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )
    end if

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
      call CloneVector ( z, x, vectorNameText='_z' )
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
    call CloneVector ( z, x, vectorNameText='_z' )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = x%quantities(i)%values - y%quantities(i)%values
    end do
  end function SubtractVectors

  ! ---------------------------------------- ValidateVectorQuantity -------
  
  ! This function performes a series of tests on a quantity to make sure it
  ! matches our requirements
  
  function ValidateVectorQuantity(quantity, coherent, stacked, regular,&
    & minorFrame, verticalCoordinate, frequencyCoordinate, &
    & noInstances, noSurfs, quantityType, molecule, sayWhyNot)

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: QUANTITY ! Test quantity
    logical, optional, intent(IN) :: COHERENT ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: STACKED  ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: REGULAR ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: MINORFRAME ! .TRUE.,.FALSE. or not present

    integer, optional, dimension(:), intent(IN) :: VERTICALCOORDINATE
    integer, optional, dimension(:), intent(IN) :: FREQUENCYCOORDINATE
    integer, optional, dimension(:), intent(IN) :: NOINSTANCES
    integer, optional, dimension(:), intent(IN) :: NOSURFS
    integer, optional, dimension(:), intent(IN) :: QUANTITYTYPE
    integer, optional, dimension(:), intent(IN) :: MOLECULE
    logical, optional, intent(IN)               :: sayWhyNot

    ! Result
    logical :: ValidateVectorQuantity

    ! Executable code

   ValidateVectorQuantity = .true.

    if (present(coherent)) then
      if (quantity%template%coherent .neqv. coherent) then
        ValidateVectorQuantity=.FALSE.
        if(present(sayWhyNot)) then
          call output('Coherent quantity checked with incoherent', advance='yes')
          call output('quantity coherent? ', advance='no')
          call output(quantity%template%coherent, advance='yes')
          call output('check coherent? ', advance='no')
          call output(coherent, advance='yes')
        end if
        return
      end if
    end if

    if (present(stacked)) then
      if (quantity%template%stacked .neqv. stacked) then
        ValidateVectorQuantity=.FALSE.
        if(present(sayWhyNot)) then
          call output('stacked quantity checked with unstacked', advance='yes')
          call output('quantity stacked? ', advance='no')
          call output(quantity%template%stacked, advance='yes')
          call output('check stacked? ', advance='no')
          call output(stacked, advance='yes')
        end if
        return
      end if
    end if

    if (present(regular)) then
      if (quantity%template%regular .neqv. regular) then
        ValidateVectorQuantity=.FALSE.
        if(present(sayWhyNot)) then
          call output('Regular quantity checked with irregular', advance='yes')
          call output('quantity regular? ', advance='no')
          call output(quantity%template%regular, advance='yes')
          call output('check regular? ', advance='no')
          call output(regular, advance='yes')
        end if
        return
      end if
    end if

    if (present(minorFrame)) then
      if (quantity%template%minorFrame .neqv. minorFrame) then
        ValidateVectorQuantity=.FALSE.
        if(present(sayWhyNot)) then
          call output('Minor fram quantity checked with not', advance='yes')
          call output('quantity minor frame? ', advance='no')
          call output(quantity%template%minorFrame, advance='yes')
          call output('check minorFrame? ', advance='no')
          call output(minorFrame, advance='yes')
        end if
        return
      end if
    end if

    if (present(verticalCoordinate)) then
      ValidateVectorQuantity=any(quantity%template%verticalCoordinate == verticalCoordinate)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with dif vert coord', advance='yes')
          call output('quantity vert coord ', advance='no')
          call output(quantity%template%verticalCoordinate, advance='yes')
          call output('check vert coord ', advance='no')
          call output(verticalCoordinate, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

    if (present(frequencyCoordinate)) then
      ValidateVectorQuantity=any(quantity%template%frequencyCoordinate == frequencyCoordinate)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with dif freq coord', advance='yes')
          call output('quantity freq coord ', advance='no')
          call output(quantity%template%frequencyCoordinate, advance='yes')
          call output('check freq coord ', advance='no')
          call output(frequencyCoordinate, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

    if (present(noInstances)) then
      ValidateVectorQuantity=any(quantity%template%noInstances == noInstances)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with dif num insts', advance='yes')
          call output('quantity num insts ', advance='no')
          call output(quantity%template%noInstances, advance='yes')
          call output('check noInstances ', advance='no')
          call output(noInstances, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

    if (present(noSurfs)) then
      ValidateVectorQuantity=any(quantity%template%noSurfs == noSurfs)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with dif num surfs', advance='yes')
          call output('quantity num surfs ', advance='no')
          call output(quantity%template%noInstances, advance='yes')
          call output('check noSurfs ', advance='no')
          call output(noSurfs, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

    if (present(quantityType)) then
      ValidateVectorQuantity=any(quantity%template%quantityType == quantityType)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with wrong type', advance='yes')
          call output('quantity type ', advance='no')
          call output(quantity%template%quantityType, advance='yes')
          call output('check quantityType ', advance='no')
          call output(quantityType, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

    if (present(molecule)) then
      ValidateVectorQuantity=any(quantity%template%molecule == molecule)
        if(present(sayWhyNot) .and. .not. ValidateVectorQuantity) then
          call output('quantity checked with wrong molecule', advance='yes')
          call output('quantity molecule ', advance='no')
          call output(quantity%template%molecule, advance='yes')
          call output('check molecule ', advance='no')
          call output(molecule, advance='yes')
        end if
      if (.not. ValidateVectorQuantity) return
    end if

  end function ValidateVectorQuantity

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
      vector%quantities(qty)%values=0.0
    end do
  end subroutine
!=======================================================================
end module VectorsModule
!=======================================================================

!
! $Log$
! Revision 2.40  2001/05/10 23:29:59  livesey
! Added some arguments to ValidateVectorQuantity
!
! Revision 2.39  2001/05/10 23:11:54  vsnyder
! Add a dumper for one vector
!
! Revision 2.38  2001/05/08 20:28:34  vsnyder
! Added stuff to dump masks
!
! Revision 2.37  2001/05/03 02:12:03  vsnyder
! Take out a line of debugging scaffolding
!
! Revision 2.36  2001/05/02 20:44:37  vsnyder
! Provide for text names for vectors that didn't come from CF
!
! Revision 2.35  2001/05/02 05:29:44  livesey
! Added index argument to GetVectorQtyByTemplateIndex
!
! Revision 2.34  2001/04/28 21:01:20  livesey
! Another bug fix in GetVectorQuantityByType
!
! Revision 2.33  2001/04/28 20:54:48  livesey
! Minor bug fix in GetVectorQuantityByType
!
! Revision 2.32  2001/04/28 07:04:32  livesey
! Minor bug fix
!
! Revision 2.31  2001/04/28 01:48:52  vsnyder
! Improve dump
!
! Revision 2.30  2001/04/28 01:27:38  livesey
! Quite an important change here.  Contents of VectorTemplate_T, VectorValue_T
! are now copies of their original entries from databases, not pointers.
!
! Revision 2.29  2001/04/25 21:57:07  livesey
! Removed insulate vector (that didn't last very long :-( )
!
! Revision 2.28  2001/04/25 20:15:23  livesey
! Tidied up InsulateVector
!
! Revision 2.27  2001/04/25 01:35:01  vsnyder
! Assignment should have been pointer assignment in CloneVector
!
! Revision 2.26  2001/04/25 01:24:54  vsnyder
! Give initial values to 'name' fields
!
! Revision 2.25  2001/04/24 21:33:53  livesey
! Added insulate vector
!
! Revision 2.24  2001/04/20 00:07:15  livesey
! Added the index field to vectorvalue_t
!
! Revision 2.23  2001/04/18 23:27:00  pwagner
! Added default .true. to Validate; also optional sayWhyNot
!
! Revision 2.22  2001/04/12 21:43:19  livesey
! Added sideband option to the quantity searches
!
! Revision 2.21  2001/04/10 22:38:20  vsnyder
! Add 'details' argument to dump routines
!
! Revision 2.20  2001/03/21 02:14:37  livesey
! Add noError argument to GetVectorQtyByType
!
! Revision 2.19  2001/03/19 17:10:47  livesey
! Added more options to validate vector quantity
!
! Revision 2.18  2001/03/16 18:17:49  livesey
! Added second vector argument and more conditions to GetVectorQuantityByType
!
! Revision 2.17  2001/03/05 00:53:59  livesey
! Added molecule argument to ValidateVectorQuantity
!
! Revision 2.16  2001/03/03 00:07:01  livesey
! Added GetVectorQtyByTemplateIndex
!
! Revision 2.15  2001/02/28 17:34:25  livesey
! Added minorFrame optional argument to ValidateVectorQuantity
!
! Revision 2.14  2001/02/27 17:18:53  livesey
! Added ValidateVectorQuantity
!
! Revision 2.13  2001/02/21 21:50:38  livesey
! Added a line to zero out a vector on creation.  Kind of like using training
! wheels in a bicycle, but avoids painful core dumps when trying to output
! unfilled vectors.
!
! Revision 2.12  2001/02/09 00:38:56  livesey
! Various changes
!
! Revision 2.11  2001/02/08 00:36:29  vsnyder
! Don't destroy in DestroyVectorValues if vector%quantities is disassociated
!
! Revision 2.10  2001/02/08 00:34:41  vsnyder
! Don't destroy in DestroyVectorInfo if vector%quantities is disassociated
!
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

