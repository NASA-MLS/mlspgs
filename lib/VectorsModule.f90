! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module VectorsModule            ! Vectors in the MLS PGS suite
!=============================================================================

  ! This module provides the simple functionality for vector quantities in the
  ! MLS Level 2 software, and related programs.

!     c o n t e n t s
!     - - - - - - - -

!         Defined types
! VectorTemplate_T   vector template (contains quantity templates)
! VectorValue_T      values corresponding to a single quantity
! Vector_T           a vector has 1 VectorTemplate_T, many VectorValue_Ts

!         Functions, operations, routines
! AddToVector                  X = X + Y
! AddVectors                   Result Z = X + Y
! AddVectorTemplateToDatabase  Adds a vector template to a database of such templates
! AddVectorToDatabase          Adds a vector to a database of such vectors
! AssignVector                 Destroy 1st arg, then assign 2nd arg to it
! AxPy                         Result z = A x + y
! ClearMask                    Clear bits of MASK according to TO_CLEAR
! ClearUnderMask               Clear elements of z corresponding to MASK
! ClearVector                  Clear elements of z
! CloneVector                  Destroy 1st arg, then use 2nd arg for a template
! ConstantXVector              Result z = A x
! ConstructVectorTemplate      Creates a vectorTemplate from a list of quantities
! CopyVector                   z = x, including copying values and mask
! CopyVectorMask               Copy mask for x to z, assuming compatible vectors
! CreateMaskArray              Allocate a MASK array
! CreateMask                   Allocate the MASK array for a vector quantity
! CreateVector                 Creates an empty vector according to a given template
! DestroyVectorDatabase        Destroys a vector database
! DestroyVectorInfo            Destroy a vector
! DestroyVectorMask            Destroys the masks stored in the vector
! DestroyVectorTemplateDatabase Destroys a vector template database
! DestroyVectorTemplateInfo    Destroys a vector template
! DestroyVectorValue           Destroy the "values" field in all of the quantities in a vector
! DotVectors                   z = x . y
! DotVectorsMasked             z = x . y, but only where mask is "off"
! DumpMask                     Display only the mask information for a vector
! dump                         Interface for next three
! dump_vector                  Display how a single vector is made up
! dump_vectors                 Display how vector database is made up
! dump_vector_templates        Display how vector template database is made up
! GetVectorQuantity            Returns pointer to quantity by name in vector
! GetVectorQuantityByType      Returns pointer to quantity by type in vector
! GetVectorQtyByTemplateIndex  Returns pointer to quantity by template in vector
! GetVectorQuantityIndexByName Returns index to quantity by name in vector
! GetVectorQuantityIndexByType Returns index to quantity by type in vector
! isVectorQtyMasked            Is the mask for VectorQty set for address
! MaskVectorQty                Set the mask for VectorQty for spec. address
! MultiplyVectors              Z = X # Y if Z present; else X = X # Y
! rmVectorFromDatabase         Removes a vector from a database of such vectors
! ScaleVector                  Y = A*X if Y is present, else X = A*X.
! SetMask                      Set bits of MASK indexed by elements of TO_SET
! SubtractFromVector           x = x - y
! SubtractVectors              Returns z = x - y
! ValidateVectorQuantity       Test vector quantity for matching components

  !---------------------------------------------------------------------------

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use DUMP_0, only: DUMP
  use Intrinsic, only: LIT_INDICES, PHYQ_INVALID
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use MLSSignals_m, only: MODULES, SIGNALS
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use STRING_TABLE, only: DISPLAY_STRING, GET_STRING, STRING_LENGTH
  use SYMBOL_TABLE, only: ENTER_TERMINAL
  use SYMBOL_TYPES, only: T_IDENTIFIER

  implicit none
  private
  ! Generics
  public :: Assignment (=), DUMP, Multiply, operator (+), operator (-)
  public :: operator (*), operator ( .DOT. ), operator ( .MDOT. )
  ! Specifics
  public :: AddToVector, AddVectors, AddVectorTemplateToDatabase
  public :: AddVectorToDatabase, AssignVector, AXPY, ClearMask
  public :: ClearUnderMask, ClearVector, CloneVector, ConstantXVector
  public :: ConstructVectorTemplate, CopyVector, CopyVectorMask, CreateMaskArray
  public :: CreateMask, CreateVector, DestroyVectorDatabase, DestroyVectorInfo
  public :: DestroyVectorMask, DestroyVectorTemplateDatabase
  public :: DestroyVectorTemplateInfo, DestroyVectorValue, DotVectors
  public :: DotVectorsMasked
  public :: DumpMask, DumpQuantityMask, DumpVectorMask, Dump_Vector
  public :: Dump_Vectors, Dump_Vector_Templates
  public :: GetVectorQuantity, GetVectorQuantityByType
  public :: GetVectorQtyByTemplateIndex, GetVectorQuantityIndexByName
  public :: GetVectorQuantityIndexByType, IsVectorQtyMasked, MultiplyVectors
  public :: RmVectorFromDatabase, ScaleVector, SetMask, SubtractFromVector
  public :: SubtractVectors, ValidateVectorQuantity, MaskVectorQty
  ! Types
  public :: VectorTemplate_T, VectorValue_T, Vector_T
  ! Parameters
  public :: M_Fill, M_FullDerivatives, M_LinAlg

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment (=)
    module procedure AssignVector
  end interface

  interface DUMP
    module procedure DUMP_VECTOR, DUMP_VECTORS, DUMP_VECTOR_TEMPLATES
  end interface

  interface DumpMask
    module procedure DumpQuantityMask, DumpVectorMask
  end interface

  interface Multiply
    module procedure MultiplyVectors
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

  interface operator ( .MDOT. )
    module procedure DotVectorsMasked
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This type describes a vector template

  type :: VectorTemplate_T
     
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
    character, dimension(:,:), pointer :: MASK => NULL() ! MASK is used to
    ! control whether elements of vectors are of interest. If MASK is not
    ! associated, every element is of interest.  Otherwise,the dimensions of
    ! MASK are the same as VALUES.  Bits of MASK(i,j) are used to determine
    ! what is not interesting.  Zero means something about VALUES(i,j) is
    ! interesting, and one means it is not.  The low-order bit is used for
    ! linear algebra.  Other bits can be used for other purposes.
  end type VectorValue_T

  ! Bit of MASK field of VectorValue_T
  integer, parameter :: M_Fill = 2**2
  integer, parameter :: M_FullDerivatives = 2**1
  integer, parameter :: M_LinAlg = 2**0

  ! This type describes a vector.

  type Vector_T
    integer :: Name = 0        ! Sub-rosa index of the vector name
    integer :: GlobalUnit = PHYQ_Invalid ! Alternative units for whole vector
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

  ! This routine adds a vector to a database of such vectors, 
  ! creating the database if necessary.

    ! Dummy arguments
    type (Vector_T), dimension(:), pointer :: DATABASE
    type (Vector_T), intent(in) ::            ITEM

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
    z%globalUnit = x%globalUnit
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
  subroutine ClearMask ( MASK, TO_CLEAR, WHAT )
  ! Clear bits of MASK indexed by elements of TO_CLEAR.  Numbering of mask
  ! elements starts at one, not zero!  If TO_CLEAR is absent, clear all of
  ! the bits of MASK.  If WHAT is absent, clear all bits.  If WHAT is
  ! present, clear only bits of MASK that correspond to "one" bits of WHAT.
    character, intent(inout), dimension(:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_CLEAR
    integer, intent(in), optional :: WHAT
    integer :: MyWhat
    MyWhat = 0
    if ( present(what) ) myWhat = not(what)
    if ( present(to_clear) ) then
      mask(to_clear) = char(iand(ichar(mask(to_clear)),myWhat))
    else
      mask = char(iand(ichar(mask),myWhat))
    end if
  end subroutine ClearMask

  !-----------------------------------------------  ClearUnderMask -----
  subroutine ClearUnderMask ( Z, Inst, Quant, What )
  ! Clear elements of Z that correspond to nonzero bits in its mask.
  ! If Inst is present, clear only elements of that instance.
  ! If Quant is present, clear only elements of that quantity.
  ! If What is present, it specifies which bits of the mask indicate which
  ! elements of Z to clear.  If What is absent, the M_LinAlg bit of the
  ! mask is used.
    type(Vector_T),intent(inout) :: Z
    integer, intent(in), optional :: Inst, Quant, What
    integer :: I1, I2    ! Bounds for instances
    integer :: II        ! Index/subscript for instances
    integer :: MyWhat    ! In case What is absent
    integer :: Q1, Q2    ! Bounds for quantities
    integer :: QI        ! Index/subscript for quantities
    integer :: VI        ! Index/subscript for values

    myWhat = m_LinAlg
    if ( present(what) ) myWhat = what
    q1 = 1
    q2 = z%template%noQuantities
    if ( present(quant) ) then
      q1 = quant
      q2 = quant
    end if
    do qi = q1, q2
      if ( associated(z%quantities(qi)%mask) ) then
        i1 = 1
        i2 = z%quantities(qi)%template%noInstances
        if ( present(inst) ) then
          i1 = inst
          i2 = inst
        end if
        do ii = i1, i2
          do vi = 1, size(z%quantities(qi)%values,1)
            if ( iand(ichar(z%quantities(qi)%mask(vi,ii)), myWhat) /= 0 ) &
              & z%quantities(qi)%values(vi,ii) = 0.0
          end do ! vi
        end do ! ii
      end if
    end do ! qi
  end subroutine ClearUnderMask 

  !-------------------------------------------------  ClearVector  -----
  subroutine ClearVector ( Z )
  ! Clear the elements of Z (make them zero).  Don't change its structure
  ! or mask.
    type(Vector_T), intent(inout) :: Z
    integer :: I
    do i = 1, size(z%quantities)
      z%quantities(i)%values = 0.0
    end do
  end subroutine ClearVector

  !-------------------------------------------------  CloneVector  -----
  subroutine CloneVector ( Z, X, VectorNameText, Database )
  ! Destroy Z, except its name.
  ! Create the characteristics of a vector to be the same template as a
  ! given one (except it has no name).  Values are allocated, but not
  ! filled.  Z's mask is allocated if X's is allocated, but it is not filled.
  ! If Database is present, add Z to it.  The position in the database
  ! isn't returned.

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to invoke DestroyVectorInfo using Z after it is no
  ! longer needed. Otherwise, a memory leak will result.  Also see
  ! AssignVector.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(Vector_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: X
    character(len=*), intent(in), optional :: VectorNameText
    type(Vector_T), dimension(:), pointer, optional :: Database
    ! Local variables:
    integer :: I, Status
    ! Executable statements:
    call destroyVectorInfo ( z )
    if ( present(vectorNameText) ) &
      & z%name = enter_terminal ( vectorNameText, t_identifier )
    z%globalUnit = x%globalUnit
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
    if ( present(database) ) i = addVectorToDatabase ( database, z )
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
  subroutine CopyVector ( Z, X, CLONE, Quant, Inst, NoValues, NoMask, &
    & VectorNameText, Database )
  ! If CLONE is present and .true., Destroy Z, deep Z = X, except the
  ! name of Z is not changed.  Otherwise, copy only the values and mask
  ! of X to Z.  If NoValues or NoMask is present and true, don't copy
  ! that part of the vector.  VectorNameText and Database are passed
  ! to CloneVector if Clone is present and true.

    type(Vector_T), intent(inout) :: Z
    type(Vector_T), intent(in) :: X
    logical, intent(in), optional :: CLONE
    integer, intent(in), optional :: Quant, Inst  ! If Quant is present,
    !  only that quantity is copied.  If furthermore Inst is present,
    !  only that instance is copied.  If Inst is present but Quant
    !  is not, the entire vector is copied.
    logical, intent(in), optional :: NoValues, NoMask  ! If present and true,
    !  don't copy the values/mask
    character(len=*), intent(in), optional :: VectorNameText
    type(Vector_T), dimension(:), pointer, optional :: Database
    logical :: DoMask, DoValues
    integer :: I
    logical MyClone
    myclone = .false.
    if ( present(clone) ) myclone = clone
    doMask = .true.
    if ( present(noMask) ) doMask = .not. noMask
    doValues = .true.
    if ( present(noValues) ) doValues = .not. noValues
    if ( myclone ) then
      call cloneVector ( Z, X, vectorNameText=vectorNameText, database=database )
    else
      if ( x%template%id /= z%template%id ) call MLSMessage &
        & ( MLSMSG_Error, ModuleName, 'Incompatible vectors in CopyVector' )
    end if
    if ( present(quant) ) then
      if ( present(inst) ) then
        if ( doValues ) z%quantities(quant)%values(:,inst) = &
            & x%quantities(quant)%values(:,inst)
        if ( doMask .and. associated (x%quantities(quant)%mask ) ) &
          z%quantities(quant)%mask(:,inst) = x%quantities(quant)%mask(:,inst)
      else
        if ( doValues ) &
          & z%quantities(quant)%values = x%quantities(quant)%values
        if ( doMask .and. associated (x%quantities(quant)%mask ) ) &
        & z%quantities(quant)%mask = x%quantities(i)%mask
      end if
    else
      do i = 1, size(x%quantities)
        if ( doValues ) z%quantities(i)%values = x%quantities(i)%values
        if ( doMask .and. associated (x%quantities(i)%mask ) ) then
          call CreateMask ( z%quantities(i) )
          z%quantities(i)%mask = x%quantities(i)%mask
        end if
      end do
    end if
  end subroutine CopyVector

  ! --------------------------------------------- CopyVectorMask ---------
  subroutine CopyVectorMask ( Z, X )
    ! Copy vector mask from x to z
    type (Vector_T), intent(inout) :: Z
    type (Vector_T), intent(in) :: X
    ! Local variables
    integer :: Q
    ! Executable code
    if ( x%template%name /= z%template%name ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, 'Incompatible vectors in CopyMask' )
    do q = 1, size(x%quantities)
      if ( associated ( x%quantities(q)%mask ) ) then
        call CreateMask ( z%quantities(q) )
        z%quantities(q)%mask = x%quantities(q)%mask
      else
        if ( associated ( z%quantities(q)%mask ) ) &
          & call Deallocate_test ( z%quantities(q)%mask, &
          & 'z%quantities(?)%mask', ModuleName )
      end if
    end do
  end subroutine CopyVectorMask

  ! ---------------------------------------------  CreateMaskArray  -----
  subroutine CreateMaskArray ( mask, values )
    ! Allocate the MASK array for a vector quantity.
    character, dimension(:,:), pointer :: MASK ! To create
    real(r8), dimension(:,:), pointer :: VALUES ! Template values
    call allocate_test ( mask, (size(values,1)), &
      & size(values,2), "MASK in CreateMaskArray", ModuleName )
    mask = char(0) ! All vector elements are interesting
  end subroutine CreateMaskArray

  ! -------------------------------------------------  CreateMask  -----
  subroutine CreateMask ( VectorValue )
  ! Allocate the MASK array for a vector quantity.
    type(VectorValue_T), intent(inout) :: VectorValue
    call createMaskArray ( vectorValue%mask, vectorValue%values )
  end subroutine CreateMask

  ! -----------------------------------------------  CreateVector  -----
  type(Vector_T) function CreateVector &
    & ( vectorName, vectorTemplate, quantities, VectorNameText, globalUnit ) &
    & result ( vector )

  ! This routine creates an empty vector according to a given template
  ! Its mask is not allocated.  Use CreateMask if one is needed.

    ! Dummy arguments
    integer, intent(in) :: vectorName   ! Sub_rosa index
    type (VectorTemplate_T), intent(in), target :: VectorTemplate ! For vector
    type (QuantityTemplate_T), dimension(:), intent(in), target :: Quantities
    character(len=*), intent(in), optional :: VectorNameText
    integer, intent(in), optional :: globalUnit

    ! Local variables
    integer :: QUANTITY                 ! Loop index
    integer :: STATUS                   ! From Allocate

    ! Executable code

    vector%name = vectorName
    if ( present(globalUnit) ) vector%globalUnit = globalUnit
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
    ! Let the destruction of the vector template take care of
    ! vector%template%quantities
    nullify ( vector%template%quantities )
    if ( .not. associated(vector%quantities) ) return
    call destroyVectorValue ( vector )
    call destroyVectorMask ( vector )
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

  !--------------------------------------------  DotVectorsMasked  -----
  real(r8) function DotVectorsMasked ( X, Y ) result (Z)
  ! Compute the inner product of two vectors.  Ignore elements masked
  ! by m_linAlg in either X or Y.

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Local variables:
    integer :: I, J, K        ! Subscripts and loop inductors
    ! Executable statements:
    if ( x%template%id /= y%template%id ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot .DOT. vectors having different templates" )
    z = 0.0_r8
    do i = 1, size(x%quantities)
      if ( associated(x%quantities(i)%mask) ) then
        if ( associated(y%quantities(i)%mask) ) then
          do j = 1, size(x%quantities(i)%values,1)
            do k = 1, size(x%quantities(i)%values,2)
              if ( iand(ior(ichar(x%quantities(i)%mask(j,k)), &
                &           ichar(y%quantities(i)%mask(j,k))), m_linAlg) == 0 ) &
                & z = z + x%quantities(i)%values(j,k) * &
                  &       y%quantities(i)%values(j,k)
            end do ! k
          end do ! j
        else
          do j = 1, size(x%quantities(i)%values,1)
            do k = 1, size(x%quantities(i)%values,2)
              if ( iand(ichar(x%quantities(i)%mask(j,k)), m_linAlg) == 0 ) &
                & z = z + x%quantities(i)%values(j,k) * &
                  &       y%quantities(i)%values(j,k)
            end do ! k
          end do ! j
        end if
      else if ( associated(y%quantities(i)%mask) ) then
          do j = 1, size(x%quantities(i)%values,1)
            do k = 1, size(x%quantities(i)%values,2)
              if ( iand(ichar(y%quantities(i)%mask(j,k)), m_linAlg) == 0 ) &
                & z = z + x%quantities(i)%values(j,k) * &
                  &       y%quantities(i)%values(j,k)
            end do ! k
          end do ! j
      else
        z = z + sum( x%quantities(i)%values * y%quantities(i)%values )
      end if
    end do
  end function DotVectorsMasked

  ! -----------------------------------------  DumpQuantityMask  -----
  subroutine DumpQuantityMask ( VectorQuantity )
    type (VectorValue_T), intent(in) :: VectorQuantity

    ! Local variables
    integer :: c                        ! Channel index
    integer :: i                        ! Instance index
    integer :: j                        ! Element index
    integer :: s                        ! Surface index
    integer :: w                        ! Line width used so far

    call output ( 'Quantity ' )
    call display_string ( vectorQuantity%template%name )

    if ( .not. vectorQuantity%template%regular ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Unable to dump mask for irregular quantities' )

    if ( .not. associated ( vectorQuantity%mask ) ) then
      call output ( ' has no mask.', advance='yes' )
    else
      call output ( '', advance='yes' )
      do i = 1, vectorQuantity%template%noInstances
        call output ( 'Instance: ' )
        call output ( i, advance='yes' )
        j = 0
        do s = 1, vectorQuantity%template%noSurfs
          call output ( 'Surface ' )
          call output ( s )
          call output ( ': ' )
          w = 13
          do c = 1, vectorQuantity%template%noChans
            if ( w > 74 ) then
              call output ( '', advance='yes' )
              call output ( '      ' )
              w = 6
            end if
            w = w + 3
            j = j + 1
            call output ( ichar(vectorQuantity%mask(j,i)), &
              & format='(z3.2)' )
          end do
          call output ( '', advance='yes' )
        end do                        ! Surface loop
      end do                          ! Instance loop
    end if                            ! Has a mask
  end subroutine DumpQuantityMask

  ! ---------------------------------------------  DumpVectorMask  -----
  subroutine DumpVectorMask ( VECTOR )
    type (Vector_T), intent(in) :: VECTOR

    ! Local variables
    integer :: q                        ! Quantity index

    ! Executable code
    call output ( 'Dumping mask for vector ' )
    call display_string ( vector%name, advance='yes' )

    do q = 1, size(vector%quantities)
      call dumpMask ( vector%quantities(q) )
    end do                              ! Loop over quantities
  end subroutine DumpVectorMask

  ! ------------------------------------------------  Dump_Vector  -----
  subroutine Dump_Vector ( VECTOR, DETAILS, NAME, &
    & QUANTITYTYPES, INSTRUMENTMODULES, SIGNAL_IDS, &
    & COHERENT, STACKED, REGULAR, MINORFRAME, MAJORFRAME, &
    & THENDITCHAFTERDUMP )

    ! dump quantities in vector according to whether they match
    ! all of the optional args: name, ..,majorframe
    ! if thenditchafterdump is present and TRUE,
    ! dump only first matching quantity
    type(Vector_T), intent(in) :: VECTOR
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump quantity values
    !                                        ! -1 Skip quantity details beyond names
    !                                        ! -2 Skip all quantity details
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    ! if the following are present, dump only quantities matching them
    integer, intent(in), optional, dimension(:)  :: QUANTITYTYPES
    integer, intent(in), optional, dimension(:)  :: INSTRUMENTMODULES
    integer, intent(in), optional, dimension(:)  :: SIGNAL_IDS
    logical, intent(in), optional                :: COHERENT
    logical, intent(in), optional                :: STACKED
    logical, intent(in), optional                :: REGULAR
    logical, intent(in), optional                :: MINORFRAME
    logical, intent(in), optional                :: MAJORFRAME
    logical, intent(in), optional                :: THENDITCHAFTERDUMP

    ! Local parameters
    integer :: J    ! Loop inductor, subscript
    integer :: MyDetails
    logical :: dumpThisQty
    logical :: myditchafterdump
    
    if ( present(thenditchafterdump) ) then
      myditchafterdump = thenditchafterdump
    else
      myditchafterdump = .false.
    endif
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
    if ( myDetails < -1 ) return
    do j = 1, size(vector%quantities)
      dumpThisQty = .true.
      if ( present (quantitytypes) ) dumpThisQty = &
        & any(vector%quantities(j)%template%quantitytype == quantitytypes)
      if ( present (instrumentmodules) ) dumpThisQty = &
        & any(vector%quantities(j)%template%instrumentmodule == instrumentmodules)
      if ( present (signal_ids) ) dumpThisQty = &
        & any(vector%quantities(j)%template%signal == signal_ids)
      if ( present (coherent) ) dumpThisQty = dumpThisQty .and. &
        & (vector%quantities(j)%template%coherent .eqv. coherent)
      if ( present (stacked) ) dumpThisQty = dumpThisQty .and. &
        & (vector%quantities(j)%template%stacked .eqv. stacked)
      if ( present (regular) ) dumpThisQty = dumpThisQty .and. &
        & (vector%quantities(j)%template%regular .eqv. regular)
      if ( present (minorFrame) ) dumpThisQty = dumpThisQty .and. &
        & (vector%quantities(j)%template%minorFrame .eqv. minorFrame)
      if ( present (majorFrame) ) dumpThisQty = dumpThisQty .and. &
        & (vector%quantities(j)%template%majorFrame .eqv. majorFrame)
      if ( dumpThisQty ) then
        call output ( j, 4 )
        call output ( "~" )
        if ( vector%quantities(j)%template%name /= 0 ) then
          call output ( ' Qty_Template_Name = ' )
          call display_string ( vector%quantities(j)%template%name )
        end if
        if ( myDetails < 0 ) then
          call output(' ', advance='yes')
          cycle
        endif
        call output ( ' noChans = ' )
        call output ( vector%quantities(j)%template%noChans, advance='no' )
        call output ( ' noSurfs = ' )
        call output ( vector%quantities(j)%template%noSurfs, advance='no' )
        call output ( ' noInstances = ' )
        call output ( vector%quantities(j)%template%noInstances, advance='no' )
        call output ( ' instanceLen = ' )
        call output ( vector%quantities(j)%template%instanceLen, advance='yes' )
        call output ( ' Qty_Template_ID = ' )
        call output ( vector%quantities(j)%template%id, advance='yes' )
        call output ( '    signal: ')
        if ( vector%quantities(j)%template%signal < 1 ) then
          call output ( '    (no database entry for this quantity) ', advance='yes')
        elseif ( signals(vector%quantities(j)%template%signal)%name < 1 ) then
          call output ( '    (no name in the database for this quantity) ', advance='yes')
        else
          call display_string ( signals(vector%quantities(j)%template%signal)%name, advance='yes' )
        endif
        call output ( '    instrumentmodule: ')
        if ( vector%quantities(j)%template%instrumentModule < 1 ) then
          call output ( '    (no database entry for this quantity) ', advance='yes')
        elseif ( vector%quantities(j)%template%instrumentModule < 1 ) then
          call output ( '    (no name in the database for this quantity) ', advance='yes')
        else
          call display_string ( modules(vector%quantities(j)%template%instrumentModule)%name, advance='yes' )
        endif
        call output ( '    (its index): ')
        call output ( vector%quantities(j)%template%instrumentmodule, advance='no')
        call output ( ' ', advance='yes')
        call output ( '  Minor Frame? (t/f): ')
        call output ( vector%quantities(j)%template%minorframe, advance='no')
        call output ( '  Major Frame? (t/f): ')
        call output ( vector%quantities(j)%template%majorframe, advance='yes')
        call output ( '  values array size is ')
        call output ( size(vector%quantities(j)%values(:,1)), advance='no')
        call output ( 'x')
        call output ( size(vector%quantities(j)%values(1,:)), advance='yes')
        if ( myDetails > 0 ) then
          call dump ( vector%quantities(j)%values, '  Elements = ' )
          if ( associated(vector%quantities(j)%mask) ) then
            call dump ( ichar(vector%quantities(j)%mask), name='Mask=', &
              & format='(z3.2)', width = 20 )
!           call dumpQuantityMask ( vector%quantities(j) )
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
        if ( myditchafterdump ) return
      end if
    end do ! j
  end subroutine Dump_Vector

  ! -----------------------------------------------  Dump_Vectors  -----
  subroutine Dump_Vectors ( VECTORS, DETAILS, NAME, &
    & QUANTITYTYPES, INSTRUMENTMODULES, SIGNAL_IDS, &
    & COHERENT, STACKED, REGULAR, MINORFRAME, MAJORFRAME, &
    & THENDITCHAFTERDUMP )
    ! dump all vectors according to whether their quantities match
    ! all of the optional args: name, ..,majorframe
    ! if thenditchafterdump is present and TRUE,
    ! dump only first matching vector
    type(Vector_T), intent(in) :: VECTORS(:)
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump quantity values
    !                                        ! -1 Skip quantity details beyond names
    !                                        ! -2 Skip all quantity details
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    ! if the following are present, dump only quantities matching them
    integer, intent(in), optional, dimension(:)  :: QUANTITYTYPES
    integer, intent(in), optional, dimension(:)  :: INSTRUMENTMODULES
    integer, intent(in), optional, dimension(:)  :: SIGNAL_IDS
    logical, intent(in), optional                :: COHERENT
    logical, intent(in), optional                :: STACKED
    logical, intent(in), optional                :: REGULAR
    logical, intent(in), optional                :: MINORFRAME
    logical, intent(in), optional                :: MAJORFRAME
    logical, intent(in), optional                :: THENDITCHAFTERDUMP

    ! Local parameters
    integer :: I
    logical :: dumpThisQty
    logical :: dumpThisVector
    integer :: J
    logical :: myditchafterdump

    if ( present(thenditchafterdump) ) then
      myditchafterdump = thenditchafterdump
    else
      myditchafterdump = .false.
    endif
    if ( size(vectors) > 1 ) then
      call output ( 'VECTORS: SIZE = ' )
      call output ( size(vectors), advance='yes' )
    end if
    do i = 1, size(vectors)
      ! Presume do not need to dump vector; hence preset to FALSE -- 
      ! becomes TRUE if wish to dump one or more quantities
      dumpThisVector = .false.
      if ( .not. associated(vectors(i)%quantities) ) then
        call output ( '(entry  ', advance='no' )
        call output ( i, advance='no' )
        call output ( '  in the vector database had been destroyed)  ', &
        & advance='yes' )
        cycle
      endif
      do j=1, size(vectors(i)%quantities)
        ! Presume need to dump quantity; hence preset to TRUE --
        ! becomes FALSE if fails to match a requirement
        dumpThisQty = .true.
        ! Check on requirements
        if ( present (quantitytypes) ) dumpThisQty = &
          & any(vectors(i)%quantities(j)%template%quantitytype == quantitytypes)
        if ( present (instrumentmodules) ) dumpThisQty = &
          & any(vectors(i)%quantities(j)%template%instrumentmodule == instrumentmodules)
        if ( present (signal_ids) ) dumpThisQty = &
          & any(vectors(i)%quantities(j)%template%signal == signal_ids)
        if ( present (coherent) ) dumpThisQty = dumpThisQty .and. &
          & (vectors(i)%quantities(j)%template%coherent .eqv. coherent)
        if ( present (stacked) ) dumpThisQty = dumpThisQty .and. &
          & (vectors(i)%quantities(j)%template%stacked .eqv. stacked)
        if ( present (regular) ) dumpThisQty = dumpThisQty .and. &
          & (vectors(i)%quantities(j)%template%regular .eqv. regular)
        if ( present (minorFrame) ) dumpThisQty = dumpThisQty .and. &
          & (vectors(i)%quantities(j)%template%minorFrame .eqv. minorFrame)
        if ( present (majorFrame) ) dumpThisQty = dumpThisQty .and. &
          & (vectors(i)%quantities(j)%template%majorFrame .eqv. majorFrame)
        dumpThisVector = dumpThisVector .or. dumpThisQty
      enddo
      if ( dumpThisVector ) then
        call output ( i, 4 )
        call output ( ': ' )
        call dump_vector ( vectors(i), details, name, &
        & quantitytypes, instrumentmodules, signal_ids, &
        & coherent, stacked, regular, minorframe, majorframe, &
        & thenditchafterdump )
        if ( myditchafterdump ) return
      endif
    end do ! i
  end subroutine Dump_Vectors

  ! --------------------------------------  Dump_Vector_Templates  -----
  subroutine Dump_Vector_Templates ( VECTOR_TEMPLATES, DETAILS )
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
  end subroutine Dump_Vector_Templates

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
      & quantityType, molecule, instrumentModule, radiometer, signal, &
      &   sideband, noError = present(otherVector) .or. myNoError)
    if ( index /= 0 ) then
      if ( present (foundInFirst) ) foundInFirst = .true.
      GetVectorQuantityByType => vector%quantities(index)
    else
      ! Can only get here if not found in first vector and noError or other
      ! vector
      if ( present (otherVector) ) then
        index = GetVectorQuantityIndexByType ( otherVector, &
          &  quantityType, molecule, instrumentModule, radiometer, signal, &
          &  sideband, noError=myNoError )
        if ( present (foundInFirst) ) foundInFirst = .false.
        if ( index /= 0 ) &
          & GetVectorQuantityByType => otherVector%quantities( index )
      end if
    end if
  end function GetVectorQuantityByType

  ! ------------------------------- GetVectorQtyByTemplateIndex --i
  function GetVectorQtyByTemplateIndex ( vector, quantityIndex, indexInVector )
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
      msg = trim(msg) // ' that has type'
      call get_string ( lit_indices(quantityType), msg(len_trim(msg)+2:) )
      call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )
    end if

  end function GetVectorQuantityIndexByType

  ! -------------------------------  IsVectorQtyMasked  -----
  logical function IsVectorQtyMasked ( vectorQty, Row, Column, What )

  ! Is the mask for VectorQty set for address (Row, Column) ?
  ! If What is present, look at the bits of the mask specified by the union
  ! of the nonzero bits of What.  Otherwise, look at the M_LinAlg bit.
  
  ! Formal args
    type (VectorValue_T), intent(in) :: vectorQty
    integer, intent(in) ::              ROW
    integer, intent(in) ::              COLUMN
    integer, intent(in), optional ::    WHAT

    integer :: MyWhat

    myWhat = m_LinAlg
    if ( present(what) ) myWhat = what
    isVectorQtyMasked = .false.
    if ( .not. associated(vectorQty%mask)) return
    isVectorQtyMasked = iand(ichar(vectorQty%mask(row, column)), myWhat) /= 0

  end function IsVectorQtyMasked

  ! -------------------------------  MaskVectorQty  -----
  subroutine MaskVectorQty ( vectorQty, Row, Column, What )

  ! Set bits of the mask for VectorQty(Row,Column); meaning
  ! If set, don't use vectorQty%values(Row, Column)
  ! Otherwise, go ahead.  If What is present, set the bits in
  ! mask indicated by What.  Otherwise, set the M_LinAlg bit.
  
  ! Formal args
    type (VectorValue_T), intent(inout) :: vectorQty
    integer, intent(in) ::              ROW
    integer, intent(in) ::              COLUMN
    integer, intent(in), optional ::    WHAT

    integer :: MyWhat

    myWhat = m_LinAlg
    if ( present(what) ) myWhat = what
    if ( .not. associated(vectorQty%mask)) return
    vectorQty%mask(row, column) = &
      & char(ior( ichar(vectorQty%mask(row, column)), myWhat ) )

  end subroutine MaskVectorQty

  !---------------------------------------------  MultiplyVectors  -----
  subroutine MultiplyVectors ( X, Y, Z, Quant, Inst )
  ! If Z is present, destroy Z and clone a new one from X, then
  ! Z = X # Y where # means "element-by-element"; otherwise X = X # Y

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    type(Vector_T), intent(in) :: Y
    type(Vector_T), intent(out), optional, target :: Z
    integer, intent(in), optional :: Quant, Inst  ! If Quant is present,
    !  only that quantity is multiplied.  If furthermore Inst is present,
    !  only that instance is multiplied.  If Inst is present but Quant
    !  is not, the entire vector is multiplied.
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
    if ( present(quant) ) then
      if ( present(inst) ) then
        result%quantities(quant)%values(:,inst) = &
          & x%quantities(quant)%values(:,inst) * &
          & y%quantities(quant)%values(:,inst)
      else
        result%quantities(quant)%values = x%quantities(quant)%values * &
          &                               y%quantities(quant)%values
      end if
    else
      do i = 1, size(x%quantities)
        result%quantities(i)%values = &
          & x%quantities(i)%values * y%quantities(i)%values
      end do
    end if
  end subroutine MultiplyVectors

  !-----------------------------------------  RmVectorFromDatabase  -----
  integer function RmVectorFromDatabase ( DATABASE, ITEM )

  ! This routine removes a vector from a database of such vectors, 
  ! deallocating the database if necessary.
  ! Alas, doesn't work--we need to know how to undecorate character tree
  ! first before we will be able to make it work; sorry (P. Wagner)

    ! Dummy arguments
    type (Vector_T), dimension(:), pointer :: DATABASE
    type (Vector_T), intent(in) ::            ITEM

    ! Local variables
    type (Vector_T), dimension(:), pointer :: tempDatabase
    logical, parameter                     :: okToDeallocEmptyDB = .FALSE.
    include "rmItemFromDatabase.f9h"
    call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot yet (ever?) rm vector from database" ) 

    rmVectorFromDatabase = newSize
  end function RmVectorFromDatabase

  !-------------------------------------------------  ScaleVector  -----
  subroutine ScaleVector ( X, A, Y )
  ! Y = A*X if Y is present, else X = A*X.

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    real(r8), intent(in) :: A
    type(Vector_T), intent(inout), optional, target :: Y
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
  subroutine SetMask ( MASK, TO_SET, MAXBIT, WHAT )
  ! Set bits of MASK indexed by elements of TO_SET.  Numbering of mask
  ! elements starts at one, not zero!  If TO_SET is absent, set bits of MASK.
  ! for all elements.  If MaxBit is present, do not set bits for any element
  ! after MaxBit.  If WHAT is present, set the bits indicated by WHAT.
  ! Otherwise, set the M_LinAlg bit.
    character, intent(inout), dimension(:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_SET
    integer, intent(in), optional :: MaxBit, What
    integer :: I, MyMaxBit, MyWhat

    myWhat = m_LinAlg
    if ( present(what) ) myWhat = what
    if ( present(to_set) ) then
      myMaxBit = huge(0)
      if ( present(maxBit) ) myMaxBit = maxBit
      do i = 1, size(to_set)
        if ( to_set(i) > 0 .and. to_set(i) <= min(myMaxBit,size(mask)) ) then
          mask(to_set(i)) = char(ior(ichar(mask(to_set(i))),myWhat))
        end if
      end do
    else
      if ( present(maxBit) ) then
        mask(:maxbit) = char(ior(ichar(mask(:maxbit)),myWhat))
      else
        mask = char(ior(ichar(mask),myWhat))
      end if
    end if
  end subroutine SetMask

  !------------------------------------------  SubtractFromVector  -----
  subroutine SubtractFromVector ( X, Y, Quant, Inst ) ! X = X - Y.

    ! Dummy arguments:
    type(Vector_T), intent(inout) :: X
    type(Vector_T), intent(in) :: Y
    integer, intent(in), optional :: Quant, Inst  ! If Quant is present,
    !  only that quantity is subtracted.  If furthermore Inst is present,
    !  only that instance is subtracted.  If Inst is present but Quant
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
  
  ! This function performs a series of tests on a quantity to make sure it
  ! matches our requirements
  
  function ValidateVectorQuantity ( quantity, coherent, stacked, regular,&
    & minorFrame, majorFrame, verticalCoordinate, frequencyCoordinate, &
    & noInstances, noSurfs, quantityType, molecule, sayWhyNot )

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: QUANTITY ! Test quantity
    logical, optional, intent(IN) :: COHERENT ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: STACKED  ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: REGULAR ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: MINORFRAME ! .TRUE.,.FALSE. or not present
    logical, optional, intent(IN) :: MAJORFRAME ! .TRUE.,.FALSE. or not present

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
          call output('Minor frame quantity checked with not', advance='yes')
          call output('quantity minor frame? ', advance='no')
          call output(quantity%template%minorFrame, advance='yes')
          call output('check minorFrame? ', advance='no')
          call output(minorFrame, advance='yes')
        end if
        return
      end if
    end if

    if (present(majorFrame)) then
      if (quantity%template%majorFrame .neqv. majorFrame) then
        ValidateVectorQuantity=.FALSE.
        if(present(sayWhyNot)) then
          call output('Major frame quantity checked with not', advance='yes')
          call output('quantity major frame? ', advance='no')
          call output(quantity%template%majorFrame, advance='yes')
          call output('check majorFrame? ', advance='no')
          call output(majorFrame, advance='yes')
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
! Revision 2.79  2002/03/08 08:06:32  livesey
! Added explicit fill mask
!
! Revision 2.78  2002/02/14 23:15:22  vsnyder
! Add .mdot. operator
!
! Revision 2.77  2002/02/08 22:58:14  livesey
! Made CopyVectorMask public
!
! Revision 2.76  2002/02/08 22:51:40  livesey
! Added CopyVectorMask
!
! Revision 2.75  2002/02/07 02:53:06  vsnyder
! Add parameter for FullDerivatives bit for mask
!
! Revision 2.74  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.73  2002/01/18 00:34:53  livesey
! Bug fix, copyVector wasn't ensuring that mask in destination was
! allocated.
!
! Revision 2.72  2001/10/23 16:39:28  pwagner
! Added MaskVectorQty
!
! Revision 2.71  2001/10/18 23:49:46  livesey
! Tidied up a floating comma in dump_vector
!
! Revision 2.70  2001/10/18 23:31:56  pwagner
! Expanded use of details in dump_vectors; stops if try to rmVectorFromDatabase
!
! Revision 2.69  2001/10/15 22:11:54  livesey
! Added globalUnit stuff
!
! Revision 2.68  2001/10/12 23:10:18  pwagner
! Better dumps, fewer bumps
!
! Revision 2.67  2001/10/09 23:43:42  pwagner
! Some further improvements in dumping vectors
!
! Revision 2.66  2001/10/08 23:40:49  pwagner
! Improved dump routines
!
! Revision 2.65  2001/10/05 17:33:55  vsnyder
! Don't set more bits in the mask than there are elements of VALUES
!
! Revision 2.64  2001/10/04 01:50:33  vsnyder
! Add 'database' argument to CloneVector, CopyVector; cosmetic changes
!
! Revision 2.63  2001/10/03 23:05:42  vsnyder
! Added 'VectorNameText' argument to CopyVector
!
! Revision 2.62  2001/10/02 23:39:15  vsnyder
! Add quantity name to error message in GetVectorQuantityIndexByType
!
! Revision 2.61  2001/10/02 19:00:50  vsnyder
! Add ClearVector subroutine
!
! Revision 2.60  2001/10/01 20:32:27  vsnyder
! Handle word and bit indexing in mask consistently
!
! Revision 2.59  2001/09/29 00:25:51  vsnyder
! Correct word indexing for mask operations
!
! Revision 2.58  2001/09/25 19:41:07  livesey
! Added DumpMask
!
! Revision 2.57  2001/09/25 00:47:08  vsnyder
! Add noMask & noValues optional arguments to CopyVector
!
! Revision 2.56  2001/09/25 00:18:23  livesey
! Bug fix
!
! Revision 2.55  2001/09/24 23:01:11  vsnyder
! Make consistent/correct lower bound calculation for MASK array
!
! Revision 2.54  2001/09/21 17:38:46  pwagner
! Added args to dump_vector(s)
!
! Revision 2.53  2001/09/20 23:02:31  vsnyder
! Specified explicitly which entities are public (so as not to re-publish
! everything gotten by USE).  Added ClearUnderMask subroutine.
!
! Revision 2.52  2001/09/20 20:56:34  pwagner
! Added contents list; tweaked some things
!
! Revision 2.51  2001/09/19 23:40:53  pwagner
! Added rmVectorFromDatabase, isVectorQtyMasked functions
!
! Revision 2.50  2001/09/17 23:10:49  pwagner
! New optional arg majorFrame in Validate..
!
! Revision 2.49  2001/07/19 17:57:15  vsnyder
! Added 'Quant' and 'Inst' arguments to CopyVector and MultiplyVectors
!
! Revision 2.48  2001/07/17 17:33:21  livesey
! Added CreateMaskArray
!
! Revision 2.47  2001/07/06 22:04:02  livesey
! Added call to DestroyVectorMask in DestroyVectorInfo
!
! Revision 2.46  2001/06/26 20:32:31  vsnyder
! Simplify mask handling by using zero origin for first dimension
!
! Revision 2.45  2001/06/01 01:04:22  vsnyder
! Add 'Multiply' generic
!
! Revision 2.44  2001/05/25 22:33:07  livesey
! Changed a comment
!
! Revision 2.43  2001/05/17 20:17:00  vsnyder
! Don't clobber Y argument of ScaleVector by making it intent(out) -- we
! need to check its template.
!
! Revision 2.42  2001/05/11 23:33:29  vsnyder
! Get rid of double-printing of 'Without mask'
!
! Revision 2.41  2001/05/11 22:01:00  vsnyder
! Simplify dumping just one vector
!
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

