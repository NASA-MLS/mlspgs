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

  interface DUMP
    module procedure DUMP_VECTORS, DUMP_VECTOR_TEMPLATES
  end interface

  interface operator (+)
    module procedure AddVectors
  end interface

  interface operator (*)
    module procedure MultiplyVectors ! element-by-element
    module procedure ScaleVector     ! by a scalar on the left
  end interface

  interface operator ( .DOT. )
    module procedure DotVectors
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
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
    integer, dimension(:), pointer :: Starts ! Where the quantities' values
      ! start in the rank-one quantity vector -- see Vector_T below.

    ! Indices of the quantity templates in the quantities database

    integer, dimension(:), pointer :: QUANTITIES
  end type VectorTemplate_T

  ! This type describes the subset of the values of a vector that
  ! correspond to a single quantity.

  type VectorValue_T
    type (QuantityTemplate_T), pointer :: TEMPLATE ! Template for this quantity
    ! The dimensions of VALUES are Frequencies (or 1), Vertical Coordinates
    ! (or 1), and Horizontal Instances (scan or profile or 1).  The target of
    ! VALUES is the appropriate piece of the VALUES field of the containing
    ! object of type Vector_T.
    real(r8), dimension(:,:,:), pointer :: VALUES => NULL()
  end type VectorValue_T

  ! This type describes a vector.

  type Vector_T
    integer :: Name            ! Sub-rosa index of the vector name
    type (VectorTemplate_T), pointer :: TEMPLATE ! In the template database
    ! There is a subterfuge going on here.  We allocate VALUES, and then
    ! take pointers to pieces of it as QUANTITIES%VALUES.  This seems to
    ! be impossible, because VALUES is rank 1, but QUANTITIES%VALUES are
    ! rank 3.  We do this by way of an external subroutine, GET_3D_VIEW,
    ! that is declared here to have a rank-1 first argument, but is declared
    ! within itself to have a rank-3 first argument.  That is WE LIE.  This
    ! is OK, because rank-1 vectors are required to be stored contiguously.
    ! This all worked in Fortran 77.  Hopefully, Fortran will someday have
    ! a straightforward way to take a rank-3 view of a rank-1 vector.
    real(r8), dimension(:), pointer :: VALUES => NULL()
    ! The dimension of QUANTITIES is the same as for the QUANTITIES field
    ! of the vector template.  Each element of QUANTITIES here corresponds
    ! to the one in the same position of the QUANTITIES field of the
    ! vector template.
    type (VectorValue_T), dimension(:), pointer :: QUANTITIES
  end type Vector_T

  ! This incrementing counter is used to set the id field for a vector template

  integer, save, private :: vectorTemplateCounter = 0

  interface 
    subroutine GET_3D_VIEW ( Input, I1, I2, I3, Output )
      use MLSCommon, only: R8
      integer, intent(in) :: I1, I2, I3
      real(r8), intent(in), target :: Input(*)
      real(r8), pointer :: Output(:,:,:)
    end subroutine GET_3D_VIEW
  end interface

contains ! =====     Public Procedures     =============================

  !--------------------------------------------------  AddVectors  -----
  type (Vector_T) function AddVectors ( X, Y ) result (Z)
  ! Add two vectors, producing one having the same template (but no name,
  ! of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call CloneVector ( z, x )
    z%values = x%values + y%values
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

    database(newSize) = item
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

    database(newSize) = item
    AddVectorToDatabase = newSize
  end function AddVectorToDatabase

  !--------------------------------------------------------  AXPY  -----
  type (Vector_T) function AXPY ( A, X, Y ) result (Z)
  ! Multiply the vector X by A and add Y to it, producing one having the
  ! same template (but no name, of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    real, intent(in) :: A
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call CloneVector ( z, x )
    z%values = a * x%values + y%values
  end function AXPY

  !-------------------------------------------------  CloneVector  -----
  subroutine CloneVector ( Z, X )
  ! Create the characteristics of a vector to be the same template as a
  ! given one (except it has no name).  Values are allocated, but not
  ! filled.

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(out) :: Z
    type(Vector_T), intent(in) :: X
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    z%name = 0
    z%template = x%template
    allocate ( z%quantities(size(x%quantities)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // "z%quantities" )
    do i = 1, size(x%quantities)
      z%quantities(i)%template => x%quantities(i)%template
    end do
    call createValues ( z )
  end subroutine  CloneVector

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
    integer :: qty, status

    ! Executable code
    vectorTemplate%name = name
    vectorTemplate%noQuantities = SIZE(selected)
    call allocate_test ( vectorTemplate%starts, vectorTemplate%noQuantities, &
      & "vectorTemplate%starts", ModuleName )
    vectorTemplate%starts(1) = 1
    do qty = 1, vectorTemplate%noQuantities - 1
      vectorTemplate%starts(qty+1) = vectorTemplate%starts(qty) + &
        & quantities(selected(qty))%noInstances * &
        & quantities(selected(qty))%instanceLen
    end do
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

  ! -----------------------------------------------  CreateVector  -----
  type(Vector_T) function CreateVector &
    & ( vectorName, vectorTemplate, quantities ) result (vector )

  ! This routine creates an empty vector according to a given template

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
  subroutine DestroyVectorInfo ( vector )

  ! This routine destroys a vector created above

    ! Dummy arguments
    type (Vector_T), intent(inout) :: VECTOR

    ! Local Variables
    integer :: I, Status

    ! Executable code

    call deallocate_test ( vector%values, "vector%values", ModuleName )
    deallocate ( vector%quantities, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_deallocate // "vector%quantities" )
  end subroutine DestroyVectorInfo

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
  subroutine DestroyVectorTemplateInfo ( vectorTemplate )

  ! This subroutine destroys a vector template created above

    ! Dummy arguments
    type (VectorTemplate_T), intent(inout) :: vectorTemplate

    ! Local variables
    integer :: QTY, STATUS

    ! Executable code

    deallocate ( vectorTemplate%quantities, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & MLSMSG_deallocate // "vectorTemplate%quantities" )

    vectorTemplate%noQuantities = 0
    vectorTemplate%totalInstances = 0
    vectorTemplate%totalElements = 0
    vectorTemplate%id = 0
    vectorTemplate%name = 0
  end subroutine DestroyVectorTemplateInfo

  !--------------------------------------------------  DotVectors  -----
  real(r8) function DotVectors ( X, Y ) result (Z)
  ! Compute the inner product of two vectors.

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, K, L
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot .DOT. vectors having different templates" )
    z = dot_product(x%values,y%values)
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
    integer :: I, J
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
  function GetVectorQuantity ( vector, quantity, quantityIsName )

  ! This function returns a pointer to the information about one quantity
  ! within a vector.

    ! Dummy arguments
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: Quantity                 ! Quantity index or name
    ! If a Quantity index, it indexes the Quantities field of Vector_T,
    ! not the quantities data base (which we don't have access to here,
    ! anyway).
    logical, intent(in), optional :: quantityIsName ! Quantity is Sub-rosa

    ! Result
    type(VectorValue_T), pointer :: GetVectorQuantity

    ! Local variables
    character(len=127) :: MSG
    integer :: Search

    ! Executable code
    if ( present(quantityIsName) ) then
      if ( quantityIsName ) then
        do search = 1, size(vector%quantities)
          if ( quantity == vector%quantities(search)%template%name ) then
            GetVectorQuantity => vector%quantities(search)
          end if
        end do
        call get_string ( quantity, msg )
        msg(string_length(quantity)+2:) = 'is not a quantity in vector'
        call get_string ( vector%name, msg(len_trim(msg)+2:) )
        call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )
      end if
    end if

    if ( quantity < 1 .or. quantity > size(vector%quantities) ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & '"Quantity" is out of range as a subscript for "vector%quantities".' )
    GetVectorQuantity => vector%quantities(quantity)

  end function GetVectorQuantity

  !---------------------------------------------  MultiplyVectors  -----
  type (Vector_T) function MultiplyVectors ( X, Y ) result (Z)
  ! Multiply two vectors element-by-element, producing one having the
  ! same template (but no name, of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot multiply vectors having different templates" )
    call CloneVector ( z, x )
    z%values = x%values * y%values
  end function MultiplyVectors

  !-------------------------------------------------  ScaleVector  -----
  type (Vector_T) function ScaleVector ( A, X ) result (Z)
  ! Multiply the vector X by A, producing one having the same template
  ! (but no name, of course).

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using the result of this
  ! function after it is no longer needed. Otherwise, a memory leak will
  ! result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    real(r8), intent(in) :: A
    type(Vector_T), intent(in) :: X
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    call CloneVector ( z, x )
    z%values = a * x%values
  end function ScaleVector

! =====     Private Procedures     =====================================
  subroutine CreateValues ( Vector )
  ! Allocate space for the values of a vector.  Create rank-3 views of
  ! the values for each of the quantities.
    type(Vector_T), intent(inout) :: Vector
    integer :: QTY
    call allocate_test ( vector%values, vector%template%totalElements, &
      & "vector%values", ModuleName )
    do qty = 1, size(vector%quantities)
      call get_3d_view ( vector%values(vector%template%starts(qty):), &
        & vector%quantities(qty)%template%noChans, &
        & vector%quantities(qty)%template%noSurfs, &
        & vector%quantities(qty)%template%noInstances, &
        & vector%quantities(qty)%values )
    end do
  end subroutine
!=======================================================================
end module VectorsModule
!=======================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

