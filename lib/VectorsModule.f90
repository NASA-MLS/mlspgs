! Copyright 2012, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module VectorsModule            ! Vectors in the MLS PGS suite
!=============================================================================

  ! This module provides the simple functionality for vector quantities in the
  ! MLS Level 2 software, and related programs.

! === (start of toc) ===                                                 
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
! AreEqual                     Check that two vectors are equal valuewise
!                               (assumes the same template, ignores masks)
! AreUnEqual                   Check that two vectors are unequal valuewise
!                               (assumes the same template, ignores masks)
! AssignVector                 Destroy 1st arg, then assign 2nd arg to it
! AxPy                         Result z = A x + y
! CheckIntegrity_Vector
! CheckIntegrity_VectorTemplate
! CheckIntegrity_VectorValue
! CheckNaN
! CheckVectorForNaN
! CheckVectorQuantityForNaN
! ClearMask                    Clear bits of MASK according to TO_CLEAR
! ClearUnderMask               Clear elements of z corresponding to MASK
! ClearVector                  Clear elements of z
! CloneVector                  Destroy 1st arg, then use 2nd arg for a template
! CloneVectorQuantity          Destroy 1st arg, then use 2nd arg for a template
! ConstantXVector              Result z = A x
! ConstructVectorTemplate      Creates a vectorTemplate from a list of quantities
! CopyVector                   z = x, including copying values and mask
! CopyVectorMask               Copy mask for x to z, assuming compatible vectors
! CreateMask                   Allocate the MASK array for a vector quantity
! CreateVector                 Creates an empty vector according to a given template
! CreateVectorValue            Creates an empty vector value according to a given template
! DestroyVectorDatabase        Destroys a vector database
! DestroyVectorInfo            Destroy a vector
! DestroyVectorMask            Destroy the masks stored in the vector
! DestroyVectorQuantityMask    Destroy the MASK stored in one vector quantity
! DestroyVectorQuantityValue   Destroy the VALUES stored in one vector quantity
! DestroyVectorTemplateDatabase Destroys a vector template database
! DestroyVectorTemplateInfo    Destroys a vector template
! DestroyVectorValue           Destroy the "values" field in all of the quantities in a vector
! Diff                         Generic for several diff_...; see the interface
! DivideVectors                Y = A / X if Y is present, else X = A / X
! DotVectors                   z = x . y
! DotVectorsMasked             z = x . y, but only where mask is "off"
! DotVectorsMaybeMasked        z = x . y, maybe masked if optional arg is true
! DumpMask                     Display only the mask information for a vector
! Dump                         Generic for several dump_...; see the interface
! DumpQuantityMask
! DumpVectorMask
! DumpVectorNorms              Dump the vector norm, or its quantity norms
! Dump_vector                  Display how a single vector is made up
! Dump_vectors                 Display how vector database is made up
! Dump_Vector_Quantity         Display a vector quantity
! Dump_vector_templates        Display how vector template database is made up
! GatherVectorQuantity         Returns a new quantity whose values are
!                                the smaller set in the hyperslab described
!                                by start, count, stride, block
! GetVectorQuantity            Returns pointer to quantity by name in vector
! GetVectorQuantityByType      Returns pointer to quantity by type in vector
! GetVectorQtyByTemplateIndex  Returns pointer to quantity by template in vector
! GetVectorQuantityIndexByName Returns index to quantity by name in vector
! GetVectorQuantityIndexByType Returns index to quantity by type in vector
! InflateVectorDatabase
! InflateVectorTemplateDatabase
! IsVectorQtyMasked            Is the mask for VectorQty set for address
! L2Norm                       L2 norm of the vector
! L2Norms                      L2 norms of the vector's quantities
! MaskVectorQty                Set the mask for VectorQty for spec. address
! MoveVectorQuantity           Move VALUES and MASK fields from one qty to another
! MultiplyVectors              Z = X # Y if Z present; else X = X # Y
! NullifyVectorTemplate
! NullifyVectorValue
! NullifyVector
! PowVector                    X = X ** Power element-by-element
! ReciprocateVector            Y = A / X if Y is present, else X = A / X -- scalar A
! RemapVectorMask              Remap MASK1 to MASK, MASK3 and MASK4
! RemapVectorValue             Remap VALUE1 to VALUES, VALUE3, and VALUE4
! ReshapeVectorValue           Reshape source values to fit destination loosely
! ReverseMask                  Reverse bits of MASK according to TO_CLEAR
! RmVectorFromDatabase         Removes a vector from a database of such vectors
! ScaleVector                  Y = A*X if Y is present, else X = A*X.
! SetMask                      Set bits of MASK indexed by elements of TO_SET
! SubtractFromVector           x = x - y
! SubtractVectors              Returns z = x - y
! ValidateVectorQuantity       Test vector quantity for matching components
! VectorMemoryInUse            Returns number of elements in all quantities
! VectorsMemoryInUse           Returns number of elements in all vectors
! === (end of toc) ===

  ! --------------------------------------------------------------------------

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, Test_Allocate, &
    & Test_Deallocate
  use HyperSlabs, only: ExtractArray
  use, Intrinsic :: Iso_C_Binding, only: C_Intptr_T, C_Loc
  use BitStuff, only: DumpBitNames, IsBitSet
  use Diff_1, only: Diff
  use Dump_0, only: Dump
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: Lit_Indices, Phyq_Invalid, L_Vmr
  use Lexer_Core, only: Where_T
  use MLSCommon, only: M_Cloud, M_Fill, M_FullDerivatives, M_Ignore, &
    & M_LinAlg, M_Spare, M_Tikhonov
  use MLSFinds, only: FindFirst, FindUnique
  use MLSKinds, only: R8, Rv
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMSG_Error, &
    & MLSMSG_Warning
  use MLSSignals_M, only: Modules, Signals, GetSignalName
  use Output_M, only: Blanks, Newline, Output
  use QuantityTemplates, only: QuantityTemplate_T, CheckIntegrity, &
    & CopyQuantityTemplate, DestroyQuantityTemplateContents, Dump, &
    & NullifyQuantityTemplate
  use String_Table, only: Display_String, Get_String_Rude=>get_String, &
    & IsStringInTable, String_Length
  use Symbol_Table, only: Enter_Terminal
  use Symbol_Types, only: T_Identifier

  implicit none
  private
  ! Generics
  public :: Assignment (=), AreEqual, AreUnEqual, Multiply
  public :: operator (+), operator (-), operator (*)
  public :: operator (==), operator (/=), operator ( .DOT. )
  public :: operator ( .MDOT. )
  public :: DIFF, DUMP
  ! Specifics
  public :: ADDTOVECTOR, ADDVECTORS, ADDVECTORTEMPLATETODATABASE
  public :: ADDVECTORTODATABASE, ASSIGNVECTOR, AXPY
  public :: CHECKINTEGRITY, CHECKNAN, CHECKVECTORFORNAN
  public :: CHECKVECTORQUANTITYFORNAN, CLEARMASK, CLEARUNDERMASK, CLEARVECTOR
  public :: CLONEVECTOR, CLONEVECTORQUANTITY, CONSTANTXVECTOR
  public :: CONSTRUCTVECTORTEMPLATE, COPYVECTOR, COPYVECTORMASK
  public :: CREATEMASK, CREATEVECTOR, CREATEVECTORVALUE, DESTROYVECTORDATABASE
  public :: DESTROYVECTORINFO, DESTROYVECTORMASK, DESTROYVECTORQUANTITYMASK
  public :: DESTROYVECTORQUANTITYVALUE, DESTROYVECTORTEMPLATEDATABASE
  public :: DESTROYVECTORTEMPLATEINFO, DESTROYVECTORVALUE
  public :: DIFFVECTORQUANTITIES
  public :: DIVIDEVECTORS
  public :: DOTVECTORS, DOTVECTORSMASKED, DOTVECTORSMAYBEMASKED
  public :: DOTVECTORQUANTITIES, DOTVECTORQUANTITIESMASKED
  public :: DOTVECTORQUANTITIESMAYBEMASKED
  public :: DUMPNICEMASKSUMMARY, DUMPMASK, DUMPQUANTITYMASK, DUMPVECTORMASK
  public :: DUMPVECTORNORMS
  public :: DUMP_VECTOR, DUMP_VECTORS, DUMP_VECTOR_QUANTITY
  public :: DUMP_VECTOR_TEMPLATE, DUMP_VECTOR_TEMPLATES
  public :: GATHERVECTORQUANTITY, GETVECTORQUANTITY, GETVECTORQUANTITYBYTYPE
  public :: GETVECTORQTYBYTEMPLATEINDEX, GETVECTORQUANTITYINDEXBYNAME
  public :: GETVECTORQUANTITYINDEXBYTYPE, INFLATEVECTORDATABASE
  public :: INFLATEVECTORTEMPLATEDATABASE, ISVECTORQTYMASKED
  public :: L2NORM, L2NORMQ, L2NORMV
  public :: MASKVECTORQTY, MOVEVECTORQUANTITY, MULTIPLYVECTORS
  public :: NULLIFYVECTORTEMPLATE, NULLIFYVECTORVALUE, NULLIFYVECTOR, POWVECTOR
  public :: QUANTITYTEMPLATE_T ! FOR FULL F95 COMPATIBILITY
  public :: RECIPROCATEVECTOR, REMAPVECTORMASK, REMAPVECTORVALUE, REVERSEMASK
  public :: RESHAPEVECTORVALUE, RMVECTORFROMDATABASE
  public :: SCALEVECTOR, SETMASK, SUBTRACTFROMVECTOR
  public :: SUBTRACTVECTORS, VALIDATEVECTORQUANTITY
  public :: VectorMemoryInUse, VectorsMemoryInUse
  ! Types
  public :: VECTORTEMPLATE_T, VECTORVALUE_T, VECTOR_T
  ! Parameters
  public :: M_IGNORE, M_CLOUD, M_FILL, M_FULLDERIVATIVES, M_LINALG, &
    & M_SPARE, M_TIKHONOV, RV

! =====     Defined Operators and Generic Identifiers     ==============

  interface Assignment (=)
    module procedure AssignVector, AssignVectorValue
  end interface

  interface AreEqual
    module procedure AreEqual_Scalar, AreEqual_Vector
  end interface

  interface AreUnEqual
    module procedure AreUnEqual_Scalar, AreUnEqual_Vector
  end interface

  interface CheckIntegrity
    module procedure CheckIntegrity_VectorValue, CheckIntegrity_VectorTemplate, &
      & CheckIntegrity_Vector
  end interface

  interface CheckNaN
    module procedure CheckVectorForNaN, CheckVectorQuantityForNaN
  end interface

  interface CLEARMASK
    module procedure CLEARMASK_1D, CLEARMASK_2D
  end interface

  interface DIFF
    module procedure DIFFVECTORQUANTITIES
  end interface

  interface DUMP
    module procedure DUMP_VECTOR, DUMP_VECTORS, DUMP_VECTOR_QUANTITY
    module procedure DUMP_VECTOR_TEMPLATE, DUMP_VECTOR_TEMPLATES
  end interface

  interface DumpMask
    module procedure DumpQuantityMask, DumpVectorMask
  end interface

  interface GETVECTORQUANTITYINDEXBYNAME
    module procedure GETVECTORQUANTITYINDEXBYNAME_CHAR
    module procedure GETVECTORQUANTITYINDEXBYNAME_SBR
  end interface

  interface L2Norm
    module procedure L2NormQ, L2NormV
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

  interface operator (==)
    module procedure AreEqual_Scalar, AreEqual_Vector
  end interface

  interface operator (/=)
    module procedure AreUnEqual_Scalar, AreUnEqual_Vector
  end interface

  interface operator ( .DOT. )
    module procedure DotVectors, DotVectorQuantities
  end interface

  interface operator ( .MDOT. )
    module procedure DotVectorsMasked, DotVectorQuantitiesMasked
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Some introductory remarks:
  ! There is a hierarchy of user-defined types
  ! * At the bottom are the VectorValue_T which correspond to specific
  !   instances of quantities like H2O or Temperature
  ! * One level up are Vector_T, built of as many VectorValue_T
  !   as you please
  ! * In parallel with these are the templates, used for holding
  !   everything each specific instance shares with the others;
  !   e.g., type, signal, geolocations
  ! * At the bottom are the QuantityTemplate_T, holding shared data
  ! * One level up are VectorTemplate_T
  ! The following diagram summarizes these relations

  !   Templates             Specific Instances
  !   ---------             ------------------
  !  VectorTemplate_T          Vector_T
  !        ^                      ^
  !        |                      |
  !  QuantityTemplate_T      VectorValue_T
  
  ! This type describes a vector template

  type :: VectorTemplate_T
     
    ! Administrative stuff
    integer :: Name = 0           ! Sub-rosa index of name, if any, else zero
    type(where_t) :: Where        ! Source_ref for creation if by L2CF

    ! General information about the vector

    integer :: NoQuantities = 0   ! Number of quantities in the vector
    integer :: TotalInstances = 0 ! Number of horizontal instances in the vector
    integer :: TotalElements = 0  ! Total of numbers of elements in the vector
    integer, dimension(:), pointer :: QUANTITIES => NULL() ! Indices of the
    !                            quantity templates in the quantities database
  end type VectorTemplate_T

  ! This type describes the subset of the values of a vector that
  ! correspond to a single quantity. 
  
  ! ---- Should have called this QuantityValue_T. ----

  type VectorValue_T
    type (QuantityTemplate_T) :: TEMPLATE ! Template for this quantity.
    integer :: Index = 0          ! Index of this quantity in the vector database
    real(rv), dimension(:), pointer, contiguous :: VALUE1 => NULL() ! The
    ! dimension of VALUE1 is Frequencies (or 1) * Vertical Coordinates or MIF
    ! (or 1) * Horizontal Instances (scan or profile or 1) * Cross-track
    ! instances.  These are taken from (template%noChans * template%noSurfs *
    ! template%noInstances).  This is the one that's allocated and deallocated.
    real(rv), dimension(:,:), pointer, contiguous :: VALUES => NULL() ! The
    ! dimensions of VALUES are Frequencies (or 1) * Vertical Coordinates or MIF
    ! (or 1), and Horizontal Instances (scan or profile or 1) * Cross-track
    ! instances.  These are taken from (template%noChans * template%noSurfs,
    ! template%noInstances).  This is a rank remapping of VALUE1.
    real(rv), dimension(:,:,:), pointer, contiguous :: VALUE3 => NULL() ! The
    ! dimensions of VALUE3 are Frequencies, Vertical Coordinates or MIF, and
    ! Horizontal Instances (MAF or profile or 1) * Cross-track instances.  These
    ! are taken from template%noChans, template%noSurfs, and
    ! template%noInstances.  This is a rank remapping of VALUE1.
    real(rv), dimension(:,:,:,:), pointer, contiguous :: VALUE4 => NULL() ! The
    ! dimensions of VALUE3 are Frequencies, Vertical Coordinates or MIF,
    ! Horizontal Instances (MAF or profile or 1), and Cross-track instances. 
    ! These are taken from template%noChans, template%noSurfs, and
    ! template%noInstances.  This is a rank remapping of VALUE1.
    character, dimension(:), pointer, contiguous :: MASK1 => NULL() ! MASK1 is
    ! the array that is allocated and deallocated.  MASK and MASK3 are rank
    ! remappings of MASK1.
    character, dimension(:,:), pointer, contiguous :: MASK => NULL() ! MASK is
    ! used to control whether elements of vectors are of interest. If MASK is
    ! not associated, every element is of interest.  Otherwise,the dimensions of
    ! MASK are the same as VALUES.  Bits of MASK(i,j) are used to determine what
    ! is not interesting.  Zero means something about VALUES(i,j) is
    ! interesting, and one means it is not.  The low-order bit is used for
    ! linear algebra.  Other bits can be used for other purposes.
    ! Actually the mask bits are gotten at from the ichar(mask)
    ! Inversely, given the integer representation of the bits, we get the mask
    ! by mask = char(int).  This is a rank remapping of MASK1.
    character, dimension(:,:,:), pointer, contiguous :: MASK3 => NULL() ! This
    ! is used for masking VALUE3, and has the same dimensions.  This is a rank
    ! remapping of MASK1.
    character, dimension(:,:,:,:), pointer, contiguous :: MASK4 => NULL() ! This
    ! is used for masking VALUE4, and has the same dimensions.  This is a rank
    ! remapping of MASK1.
    integer :: Label = 0        ! An optional label for this to be used as for
    ! example a swath name.  Often used in conjunction with the 'batch'
    ! approach to direct writes.
    character(len=40) :: AllocationName = 'None'
    ! These fields apply only to products of a neural network model
    integer, pointer, dimension(:)    :: BinNumber => NULL()
    integer, pointer, dimension(:)    :: MAF => NULL()
  end type VectorValue_T

  character(len=16), dimension(7), parameter :: maskBitNames = (/ &
    & 'linear algebra  ', 'full derivatives', 'fill            ' , &
    & 'Tikhonov        ', 'cloud           ', 'ignore          ' , &
    & 'spare           ' /)

  ! This type describes a vector.

  type Vector_T
    integer :: Name = 0        ! Sub-rosa index of the vector name
    type(where_t) :: Where     ! Source_ref for creation if by L2CF
    integer :: GlobalUnit = PHYQ_Invalid ! Alternative units for whole vector
    type (VectorTemplate_T) :: TEMPLATE ! Template for this vector
    type (VectorValue_T), dimension(:), pointer :: QUANTITIES => NULL() ! The
    ! dimension of QUANTITIES is the same as for the QUANTITIES field of the
    ! vector template.  Each element of QUANTITIES here corresponds to the
    ! one in the same position in the QUANTITIES field of the vector template.
!   contains
!     final :: DestroyVectorInfo
  end type Vector_T

  private :: CreateValues

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  AddToVector  -----
  subroutine AddToVector ( X, Y, Scale )  ! X = X + [Scale*] Y
    ! Dummy arguments:
    type(Vector_T), intent(inout) :: X
    type(Vector_T), intent(in) :: Y
    real(r8), intent(in), optional :: Scale
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    ! Executable statements:
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    if ( present(scale) ) then
      do i = 1, size(x%quantities)
        x%quantities(i)%values = x%quantities(i)%values + &
          & scale * y%quantities(i)%values
      end do
    else
      do i = 1, size(x%quantities)
        x%quantities(i)%values = x%quantities(i)%values + y%quantities(i)%values
      end do
    end if
  end subroutine AddToVector

  ! -------------------------------------------------  AddVectors  -----
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
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call nullifyVector ( z ) ! for Sun's still useless compiler
    call CloneVector ( z, x, vectorNameText='_z' )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = x%quantities(i)%values + y%quantities(i)%values
    end do
  end function AddVectors

  ! --------------------------------  AddVectorTemplateToDatabase  -----
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

  ! ----------------------------------------  AddVectorToDatabase  -----
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

  ! --------------------------------------------  AreEqual_Scalar  -----
  logical function AreEqual_Scalar ( A, C ) result ( Equal )

  ! This routine checks that all associated values in vector A are equal to
  ! scalar C.

    ! Dummy arguments
    type (Vector_T), intent(in)           ::  A
    real(rv), intent(in)                  ::  C

    ! Local variables
    type (VectorValue_T), pointer :: AQ ! a quantity
    integer :: qIndex
    ! Executable
    equal = .true. ! In case there are no quantities
    do qIndex = 1, size ( a%quantities ) 
      aq => a%quantities(qIndex)
      if ( .not. associated(aq) ) cycle
      equal = all(aq%values == c)
      if ( .not. equal ) exit
    end do

  end function AreEqual_Scalar

  ! --------------------------------------------  AreEqual_Vector  -----
  logical function AreEqual_Vector ( A, B ) result ( Equal )

  ! This routine checks that all associated values in vector A are equal
  ! to corresponding values in vector B.

    ! Dummy arguments
    type (Vector_T), intent(in)           ::  A
    type (Vector_T), intent(in)           ::  B

    ! Local variables
    type (VectorValue_T), pointer :: AQ ! a quantity
    type (VectorValue_T), pointer :: BQ ! a quantity
    integer :: qIndex
    ! Executable
    equal = .true. ! In case there are no quantities

    do qIndex = 1, size ( a%quantities ) 
      aq => a%quantities(qIndex)
      bq => b%quantities(qIndex)
      equal = associated(aq) .eqv. associated(bq)
      if ( .not. equal ) exit
      if ( .not. associated(aq) .or. .not. associated(bq) ) cycle
      equal = all(aq%values == bq%values)
      if ( .not. equal ) exit
    end do

  end function AreEqual_Vector

  ! ------------------------------------------  AreUnEqual_Scalar  -----
  logical function AreUnEqual_Scalar ( A, C ) result ( UnEqual )
    ! Dummy arguments
    type (Vector_T), intent(in)           ::  A
    real(rv), intent(in)                  ::  C
    unEqual = .not. ( A == C )
  end function AreUnEqual_Scalar

  ! ------------------------------------------  AreUnEqual_Vector  -----
  logical function AreUnEqual_Vector ( A, B ) result ( UnEqual )
    ! Dummy arguments
    type (Vector_T), intent(in)           ::  A
    type (Vector_T), intent(in)           ::  B
    unEqual = .not. ( A == B )
  end function AreUnEqual_Vector

  ! -----------------------------------------------  AssignVector  -----
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

  ! ------------------------------------------  AssignVectorValue  -----
  subroutine AssignVectorValue ( Z, X )
  ! Turn off assignment for vector quantities
    use MLSMEssageModule, only: MLSMSG_Crash, MLSMessage
    type(VectorValue_T), intent(inout) :: Z
    type(VectorValue_T), intent(in) :: X
    call MLSMessage ( MLSMSG_Crash, moduleName, "Assigning VectorValue_T is a no-no" )
    z%index = x%index ! to avoid sniveling about unused dummy arguments
  end subroutine AssignVectorValue

  ! -------------------------------------------------------  AXPY  -----
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
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot add vectors having different templates" )
    call nullifyVector ( z ) ! for Sun's still useless compiler
    call CloneVector ( z, x, vectorNameText='_z' )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = &
        & a * x%quantities(i)%values + y%quantities(i)%values
    end do
  end function AXPY

  ! ------------------------------------------ CheckIntegrity_Vector -----------
  logical function CheckIntegrity_Vector ( vector, noError )
    type ( Vector_T ), intent(in) :: VECTOR
    logical, optional, intent(in) :: NOERROR

    ! Local variables
    integer :: MESSAGETYPE
    character (len=132) :: NAME
    integer :: QTY                      ! Loop counter
    integer :: totalInstances
    integer :: totalElements

    ! Executable code
    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    if ( vector%name /= 0 ) then
      call get_string ( vector%name, name, strip=.true. )
    else
      name = '<no name>'
    end if

    call output ( 'Checking integrity for vector '//trim(name), advance='yes' )
    CheckIntegrity_Vector = CheckIntegrity ( vector%template, .true. )

    if ( .not. associated ( vector%quantities ) .and. &
      & vector%template%noQuantities > 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector '//trim(name)//' does not have quantities associated' )
      CheckIntegrity_Vector = .false.
    end if

    totalInstances = 0
    totalElements = 0
    do qty = 1, vector%template%noQuantities
      call output ( 'Checking integrity for quantity ' )
      call display_string ( vector%quantities(qty)%template%name, &
        & strip=.true., advance='yes' )
      CheckIntegrity_Vector = CheckIntegrity_Vector .and. &
        & CheckIntegrity ( vector%quantities(qty), .true. )
      totalInstances = totalInstances + vector%quantities(qty)%template%noInstances
      totalElements = totalElements + vector%quantities(qty)%template%noInstances * &
        & vector%quantities(qty)%template%instanceLen
    end do

    if ( totalElements /= vector%template%totalElements ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector '//trim(name)//' does not have the right totalElements' )
      CheckIntegrity_Vector = .false.
    end if
    if ( totalInstances /= vector%template%totalInstances ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector '//trim(name)//' does not have the right totalInstances' )
      CheckIntegrity_Vector = .false.
    end if

  end function CheckIntegrity_Vector

  ! ------------------------------------------ CheckIntegrity_VectorTemplate ---
  logical function CheckIntegrity_VectorTemplate ( template, noError )
    type ( VectorTemplate_T), intent(in) :: TEMPLATE
    logical, optional, intent(in) :: NOERROR

    ! Local variables
    integer :: MESSAGETYPE
    character (len=132) :: NAME

    ! Executable code
    CheckIntegrity_VectorTemplate = .true.

    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    if ( template%name /= 0 ) then
      call get_string ( template%name, name, strip=.true. )
    else
      name = '<no name>'
    end if

    ! Check stuff
    if ( .not. associated ( template%quantities ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector template '//trim(name)//' does not have quantities associated' )
      CheckIntegrity_VectorTemplate = .false.
    end if

    if ( lbound ( template%quantities, 1 ) /= 1 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector template '//trim(name)//' has the wrong lbound for quantities' )
      CheckIntegrity_VectorTemplate = .false.
    end if

    if ( ubound ( template%quantities, 1 ) /= template%noQuantities ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector template '//trim(name)//' has the wrong ubound for quantities' )
      CheckIntegrity_VectorTemplate = .false.
    end if

    if ( any ( template%quantities < 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector template '//trim(name)//' contains a bad quantity index' )
      CheckIntegrity_VectorTemplate = .false.
    end if

  end function CheckIntegrity_VectorTemplate

  ! ------------------------------------------ CheckIntegrity_VectorValue ---
  logical function CheckIntegrity_VectorValue ( value, noError ) 
    type ( VectorValue_T), intent(in) :: VALUE
    logical, optional, intent(in) :: NOERROR

    ! Local variables
    integer :: MESSAGETYPE
    character (len=132) :: NAME

    ! Executable code
    ! Executable code
    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    ! Check the integrity of the template
    CheckIntegrity_VectorValue = CheckIntegrity ( value%template, .true. )

    if ( value%template%name /= 0 ) then
      call get_string ( value%template%name, name, strip=.true. )
    else
      name = '<no name>'
    end if

    ! Check the integrity of the values
    if ( .not. associated ( value%values ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector_value for '//trim(name)//' does not have <values> associated' )
      CheckIntegrity_VectorValue = .false.
    end if

    if ( any ( lbound ( value%values ) /= (/1,1/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector_value for '//trim(name)//' has the wrong lbound for values' )
      CheckIntegrity_VectorValue = .false.
    end if
    if ( any ( ubound ( value%values ) /= &
      & (/value%template%instanceLen,value%template%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The vector_value for '//trim(name)//' has the wrong ubound for values' )
      CheckIntegrity_VectorValue = .false.
    end if

    if ( associated ( value%mask ) ) then
      if ( any ( lbound ( value%mask ) /= (/1,1/) ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'The vector_value for '//trim(name)//' has the wrong lbound for mask' )
        CheckIntegrity_VectorValue = .false.
      end if
      if ( any ( ubound ( value%mask ) /= &
        & (/value%template%instanceLen,value%template%noInstances/) ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'The vector_value for '//trim(name)//' has the wrong ubound for mask' )
        CheckIntegrity_VectorValue = .false.
      end if
    end if

  end function CheckIntegrity_VectorValue

  ! ------------------------------------------  CheckVectorForNaN  -----
  logical function CheckVectorForNaN ( Vector, Print, Name )
  ! Check whether a vector has any NaNs in any of its VALUES, returning
  ! TRUE if so.
  ! This doesn't check the quantity templates.
    type (Vector_t), intent(in) :: Vector
    integer, intent(in) :: Print ! <= 0 => No printing
                                 !  > 0 => Call dump with details = print - 1
                                 !         on offending vector quantities
    character(len=*), intent(in), optional :: Name

    integer :: I

    CheckVectorForNaN = .false.
    if ( .not. associated(vector%quantities) ) return
    do i = 1, size(vector%quantities)
      CheckVectorForNaN = CheckVectorForNaN .or. &
        & checkNaN(vector%quantities(i), print,name,vector)
    end do
  end function CheckVectorForNaN

  ! -----------------------------------  CheckVectorQuantityForNaN  -----
  logical function CheckVectorQuantityForNaN ( VectorQuantity, Print, Name, Vector )
  ! Check whether a vector quantity has any NaNs in any of its VALUES, returning
  ! TRUE if so.
  ! This doesn't check the quantity templates.
    use IEEE_Arithmetic, only: IEEE_Is_Nan
    type (VectorValue_t), intent(in) :: VectorQuantity
    integer, intent(in) :: Print ! <= 0 => No printing
                                 !  > 0 => Call dump with details = print - 1
    character(len=*), intent(in), optional :: Name
    type (Vector_t), intent(in), optional :: Vector ! to get its name in case of dump

    CheckVectorQuantityForNaN = .false.
    if ( .not. associated(vectorQuantity%values) ) return
    CheckVectorQuantityForNaN = any(ieee_is_nan(vectorQuantity%values))
    if ( CheckVectorQuantityForNaN ) then
      if ( print > 0 ) &
        & call dump ( vectorQuantity, details=print-1, name=name, vector=vector )
    end if
  end function CheckVectorQuantityForNaN

  ! --------------------------------------------------  ClearMask  -----
  subroutine ClearMask_1d ( MASK, TO_CLEAR, WHAT )
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
  end subroutine ClearMask_1d

  subroutine ClearMask_2d ( MASK, TO_CLEAR, WHAT )
  ! Clear bits of MASK indexed by elements of TO_CLEAR.  Numbering of mask
  ! elements starts at one, not zero!  If TO_CLEAR is absent, clear all of
  ! the bits of MASK.  If WHAT is absent, clear all bits.  If WHAT is
  ! present, clear only bits of MASK that correspond to "one" bits of WHAT.
    character, intent(inout), dimension(:,:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_CLEAR
    integer, intent(in), optional :: WHAT
    integer :: instance
    do instance=1, size(mask, 2)
      call ClearMask ( mask(:, instance), to_clear, what )
    end do
  end subroutine ClearMask_2d

  ! ----------------------------------------------  ClearUnderMask -----
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
              & z%quantities(qi)%values(vi,ii) = 0.0_rv
          end do ! vi
        end do ! ii
      end if
    end do ! qi
  end subroutine ClearUnderMask 

  ! ------------------------------------------------  ClearVector  -----
  subroutine ClearVector ( Z, value )
  ! Clear the elements of Z (make them zero).  Don't change its structure
  ! or mask.
    type(Vector_T), intent(inout) :: Z
    real(rv), optional, intent(in) :: VALUE ! Value to set in (default zero)
    real(rv) :: MYVALUE
    integer :: I
    myValue = 0.0_rv
    if ( present ( value ) ) myValue = value
    do i = 1, size(z%quantities)
      z%quantities(i)%values = myValue
    end do
  end subroutine ClearVector

  ! ------------------------------------------------  CloneVector  -----
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
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, Status
    ! Executable statements:
    call destroyVectorInfo ( z )
    if ( present(vectorNameText) ) &
      & z%name = enter_terminal ( vectorNameText, t_identifier )
    z%globalUnit = x%globalUnit
    z%template = x%template
    allocate ( z%quantities(size(x%quantities)), stat=status )
    addr = 0
    if ( status == 0 ) then
      if ( size(z%quantities) > 0 ) addr = transfer(c_loc(z%quantities(1)), addr)
    end if
    call test_allocate ( status, moduleName, "z%quantities", &
      & uBounds=[size(x%quantities)], elementSize=storage_size(z%quantities) / 8, &
      & address=addr )
    z%quantities%index = x%quantities%index
    do i = 1, size(x%quantities)
      z%quantities(i)%template = x%quantities(i)%template
    end do
    call createValues ( z )
    do i = 1, size(x%quantities)
      if ( associated(x%quantities(i)%mask) ) then
        call createMask ( z%quantities(i) )
      end if
    end do
    if ( present(database) ) i = addVectorToDatabase ( database, z )
  end subroutine  CloneVector

  ! ------------------------------------------------  CloneVectorQuantity  -----
  subroutine CloneVectorQuantity ( Z, X, OPTIONS )
  ! Create a vector quantity to be the same template as a
  ! given one.  If the original's values are allocated, allocate the clone's
  ! and fill with the original's. Same with mask.
  
  ! If options is present and contains the string 'd', do a deep copy
  ! of the template. Otherwise, content yourself with a shallow one

  ! !!!!! ===== IMPORTANT NOTE ===== !!!!!
  ! It is important to deallocate Z's values and mask after they are no
  ! longer needed. Otherwise, a memory leak will result.
  ! !!!!! ===== END NOTE ===== !!!!! 

    ! Dummy arguments:
    type(vectorValue_T), intent(inout) :: Z
    type(vectorValue_T), intent(in) :: X
    character(len=*), optional, intent(in) :: OPTIONS
    ! Internal variables
    character(len=8) :: MYOPTIONS
    ! Executable statements:
    myOptions = ' '
    if ( present(options) ) myOptions = options
    call NullifyVectorValue ( z )
    if ( index( myOptions, 'd' ) > 0 ) then
      call CopyQuantityTemplate ( z%template, x%template )
    else
      z%template = x%template
    end if
    z%index = x%index
    if ( associated(x%values) ) then
      call createVectorValue ( z, 'z%values' )
      z%values = x%values
    end if
    if ( associated(x%mask) ) then
      call createMask ( z )
      z%mask = x%mask
    end if
    if ( associated(x%BinNumber) ) then
      allocate ( z%BinNumber(x%template%NoInstances) )
      z%BinNumber = x%BinNumber
    end if
    if ( associated(x%MAF) ) then
      allocate ( z%MAF(x%template%NoInstances) )
      z%MAF = x%MAF
    end if
    z%label = x%label
  end subroutine CloneVectorQuantity

  ! --------------------------------------------  ConstantXVector  -----
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
    call nullifyVector ( z ) ! for Sun's still useless compiler
    call CloneVector ( z, x, vectorNameText='_z' )
    do i = 1, size(x%quantities)
      z%quantities(i)%values = a * x%quantities(i)%values
    end do
  end function ConstantXVector

  ! ------------------------------------  ConstructVectorTemplate  -----
  subroutine ConstructVectorTemplate ( name, quantities, selected, &
    & vectorTemplate, where, ForWhom )

  ! This subroutine creates a vectorTemplate from a list of quantities.
  ! The default ordering is currently by quantity.  Later versions may
  ! have optional parameters to request other orderings.

    ! Dummy arguments
    integer, intent(in) :: NAME         ! Sub-rosa of vector template name
    type (QuantityTemplate_T), intent(in) :: quantities(:)
    integer, intent(in) :: selected(:)  ! Which quantities are selected?
    type (VectorTemplate_T), intent(out) :: vectorTemplate
    type(where_t), intent(in), optional :: where ! source_ref if created by L2CF
    character(len=*), intent(in), optional :: ForWhom ! use instead of ModuleName

    ! Executable code
!   Don't do the following.  The caller puts the actual argument into a database
!   using a shallow copy.  Destroying it would clobber a database item.
!   call destroyVectorTemplateInfo ( vectorTemplate ) ! avoid memory leaks

    vectorTemplate%name = name
    vectorTemplate%noQuantities = SIZE(selected)
    vectorTemplate%totalInstances = SUM(quantities(selected)%noInstances)
    vectorTemplate%totalElements = &
      & SUM(quantities(selected)%noInstances*quantities(selected)%instanceLen)
    if ( present(where) ) vectorTemplate%where = where
    
    ! Allocate some arrays
    if ( present(forWhom) ) then
      call allocate_test ( vectorTemplate%quantities, vectorTemplate%noQuantities, &
        & 'vectorTemplate%quantities', forWhom )
    else
      call allocate_test ( vectorTemplate%quantities, vectorTemplate%noQuantities, &
        & 'vectorTemplate%quantities', moduleName )
    end if

    ! Copy quantities over
    vectorTemplate%quantities = selected

  end subroutine ConstructVectorTemplate

  ! -------------------------------------------------  CopyVector  -----
  subroutine CopyVector ( Z, X, CLONE, Quant, Inst, NoValues, NoMask, &
    & VectorNameText, Database, AllowNameMismatch )
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
    logical, intent(in), optional :: AllowNameMismatch ! Trust that vectors are compatible
    logical :: DoMask, DoValues
    integer :: I
    logical MyClone
    logical MyAllow
    myclone = .false.
    if ( present(clone) ) myclone = clone
    doMask = .true.
    if ( present(noMask) ) doMask = .not. noMask
    doValues = .true.
    if ( present(noValues) ) doValues = .not. noValues
    myAllow = .false.
    if ( present(allowNameMismatch) ) myAllow = allowNameMismatch
    if ( myclone ) then
      call cloneVector ( Z, X, vectorNameText=vectorNameText, database=database )
    else
      if ( x%template%name /= z%template%name .and. .not. myAllow ) call MLSMessage &
        & ( MLSMSG_Error, ModuleName, 'Incompatible vectors in CopyVector' )
    end if
    if ( present(quant) ) then
      if ( doMask .and. .not. associated ( z%quantities(quant)%mask ) ) &
        & call CreateMask ( z%quantities(quant) )
      if ( present(inst) ) then
        if ( doValues ) z%quantities(quant)%values(:,inst) = &
            & x%quantities(quant)%values(:,inst)
        if ( doMask .and. associated (x%quantities(quant)%mask ) ) &
          & z%quantities(quant)%mask(:,inst) = x%quantities(quant)%mask(:,inst)
      else
        if ( doValues ) &
          & z%quantities(quant)%values = x%quantities(quant)%values
        if ( doMask .and. associated (x%quantities(quant)%mask ) ) &
          & z%quantities(quant)%mask = x%quantities(i)%mask
      end if
    else
      do i = 1, size(x%quantities)
        z%quantities(i)%index = i
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
        call destroyVectorQuantityMask ( z%quantities(q) )
      end if
    end do
  end subroutine CopyVectorMask

  ! -------------------------------------------------  CreateMask  -----
  subroutine CreateMask ( Value, ForWhom )
  ! Allocate the MASK array for a vector quantity.
    type(VectorValue_T), intent(inout) :: Value
    character(len=*), intent(in), optional :: ForWhom
    if ( present(forWhom) ) then
      call allocate_test ( value%mask1, size(value%values), &
        & trim(forWhom) // "%MASK1", ModuleName, fill=char(0) )
    else
      call allocate_test ( value%mask1, size(value%values), "MASK1", &
        & ModuleName, fill=char(0) )
    end if
    call remapVectorMask ( value )
  end subroutine CreateMask

  ! -----------------------------------------------  CreateVector  -----
  type(Vector_T) function CreateVector &
    & ( vectorName, vectorTemplate, quantities, VectorNameText, globalUnit, &
    & highBound, lowBound, noValues, where ) &
    & result ( vector )

    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

  ! This routine creates an empty vector according to a given template
  ! Its mask is not allocated.  Use CreateMask if one is needed.

    ! Dummy arguments
    integer, intent(in) :: vectorName   ! Sub_rosa index
    type (VectorTemplate_T), intent(in), target :: VectorTemplate ! For vector
    type (QuantityTemplate_T), dimension(:), intent(in), target :: Quantities
    character(len=*), intent(in), optional :: VectorNameText
    integer, intent(in), optional :: globalUnit
    logical, intent(in), optional :: highBound
    logical, intent(in), optional :: lowBound
    logical, intent(in), optional :: noValues ! Don't create values for it.
    type(where_t), intent(in), optional :: where    ! source_ref

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: QUANTITY                 ! Loop index
    integer :: Me = -1                  ! String index for tracing
    integer :: STATUS                   ! From Allocate
    logical :: MYNOVALUES               ! Copy of novalues

    ! Executable code

    call trace_begin ( me, "CreateVector", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

    vector%name = vectorName
    if ( present(globalUnit) ) vector%globalUnit = globalUnit
    if ( present(vectorNameText) ) &
      & vector%name = enter_terminal ( vectorNameText, t_identifier )
    vector%template = vectorTemplate
    allocate ( vector%quantities(vectorTemplate%noQuantities), STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( vectorTemplate%noQuantities > 0 ) addr = transfer(c_loc( &
        & vector%quantities(1)), addr)
    end if
    call test_allocate ( status, moduleName, "Vector quantities", &
      & uBounds=[vectorTemplate%noQuantities], &
      & elementSize = storage_size(vector%quantities) / 8, address=addr )
    do quantity = 1, vectorTemplate%noQuantities
      vector%quantities(quantity)%index = quantity
      vector%quantities(quantity)%template = &
        & quantities(vectorTemplate%quantities(quantity))
    end do
    myNoValues = .false.
    if ( present ( noValues ) ) myNoValues = noValues
    if ( .not. myNoValues) &
      & call createValues ( vector, highBound, lowBound )
    if ( present(where) ) vector%where = where

    call trace_end ( "CreateVector", cond=toggle(gen) .and. levels(gen) > 2 )

  end function CreateVector

  ! ------------------------------------------  CreateVectorValue  -----
  subroutine CreateVectorValue ( Value, What, Where )
  ! Allocate the Value1 component of Value and remap it.
  ! The Mask* components are destroyed.  Use CreateMask if you need to.
    type(vectorValue_t) :: Value
    character(len=*), intent(in) :: What
    character(len=*), intent(in), optional :: Where
    ! Local variables
    integer :: valueSize
    ! Executable
    ! In case it is already allocated, must deallocate Value
    ! before changing its allocation name
    call DestroyVectorQuantityValue ( value, &
      & forWhom=trim(moduleName) // '%CreateVectorValue' )
    valueSize = max( 1, value%template%noChans * &
        & value%template%noSurfs * &
        & value%template%noInstances * &
        & value%template%noCrossTrack )
    call destroyVectorQuantityMask ( value )
    if ( present(where) ) then
      value%AllocationName = trim(what) // "%values"
    else
      value%AllocationName = 'VALUE1'
    end if
    call allocate_test ( value%value1, &
      & valueSize, &
      & trim(value%AllocationName), moduleName )
    call remapVectorValue ( value )
  end subroutine CreateVectorValue

  ! --------------------------------------  DestroyVectorDatabase  -----
  subroutine DestroyVectorDatabase ( database )

  ! This subroutine destroys a vector database

    ! Dummy argument
    type (Vector_T),  dimension(:), pointer :: database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: l2gpIndex, S, Status

    if ( associated(database) ) then
      do l2gpIndex = 1, SIZE(database)
        call DestroyVectorInfo(database(l2gpIndex))
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, moduleName, 'database', s, address=addr )
    end if
    allocate ( database(0), stat=status )
    call test_allocate ( status, moduleName, "database" )
  end subroutine DestroyVectorDatabase

  ! ------------------------------------------  DestroyVectorInfo  -----
  subroutine DestroyVectorInfo ( Vector )

  ! This routine destroys a vector created above

    ! Dummy arguments
    type (Vector_T), intent(inout) :: VECTOR

    ! Local Variables:
    integer(c_intptr_t) :: Addr         ! For tracing
    character(len=80)   :: nameStr
    integer :: S, STATUS
    logical, parameter :: deeBug = .false.
    ! Executable code
    nameStr = ' '
    if ( vector%name > 0 ) &
      & call get_string ( vector%name, nameStr )
    if ( deeBug .and. len_trim(nameStr) > 0 ) &
      & call output ( 'Destroying Vector ' // trim(nameStr), advance='yes' )
    vector%name = 0
    vector%where = where_t(0,0)

    if ( .not. associated(vector%quantities) ) return
    if ( deeBug ) call output ( 'Destroying VectorValue', advance='yes' )
    call destroyVectorValue ( vector )
    if ( deeBug ) call output ( 'Destroying VectorMask', advance='yes' )
    call destroyVectorMask ( vector )
    s = size(vector%quantities) * storage_size(vector%quantities) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(vector%quantities(1)), addr)
    if ( deeBug ) call output ( 'Deallocating its quantities', advance='yes' )
    deallocate ( vector%quantities, stat=status )
    call test_deallocate ( status, moduleName, 'vector%quantities', s, address=addr )
  end subroutine DestroyVectorInfo

  ! ------------------------------------------  DestroyVectorMask  -----
  subroutine DestroyVectorMask ( Vector )

  ! This routine destroys the masks stored in the vector.

    ! Dummy arguments
    type (Vector_T), intent(inout) :: VECTOR

    ! Local Variables:
    integer :: I
    do i = 1, size(vector%quantities)
      call destroyVectorQuantityMask ( vector%quantities(i) )
    end do
  end subroutine DestroyVectorMask

  ! ----------------------------------  DestroyVectorQuantityMask  -----
  subroutine DestroyVectorQuantityMask ( Value, ForWhom )

    ! This routine destroys the MASK stored in one vector quantity.

    ! Dummy arguments
    type (vectorValue_t), intent(inout) :: Value
    character(len=*), intent(in), optional :: ForWhom
    if ( present(forWhom) ) then
      call deallocate_test ( value%mask1, trim(forWhom) // "%MASK1", moduleName )
    else
      call deallocate_test ( value%mask1, "MASK1", moduleName )
    end if
    nullify ( value%mask, value%mask3 )
  end subroutine DestroyVectorQuantityMask

  ! ---------------------------------  DestroyVectorQuantityValue  -----
  subroutine DestroyVectorQuantityValue ( VALUE, &
    & DESTROYMASK, DESTROYTEMPLATE, FORWHOM )

    ! This routine destroys the VALUES stored in one vector quantity.

    ! Dummy arguments
    type (vectorValue_t), intent(inout)    :: VALUE
    logical, intent(in), optional          :: DESTROYMASK
    logical, intent(in), optional          :: DESTROYTEMPLATE
    character(len=*), intent(in), optional :: FORWHOM
    ! Executable
    call deallocate_test ( value%value1, trim(value%allocationName), moduleName )
    nullify ( value%values, value%value3 )
    if ( present(destroyMask) ) then
      if ( destroyMask ) call destroyVectorQuantityMask ( value, forWhom )
    end if
    if ( present(destroyTemplate) ) then
      if ( destroyTemplate ) call DestroyQuantityTemplateContents ( value%template )
    end if
    if ( associated(value%BinNumber) ) then
      deallocate(value%BinNumber)
      nullify(value%BinNumber)
    endif
    if ( associated(value%MAF) ) then
      deallocate(value%MAF)
      nullify(value%MAF)
    endif
  end subroutine DestroyVectorQuantityValue

  ! ------------------------------  DestroyVectorTemplateDatabase  -----
  subroutine DestroyVectorTemplateDatabase ( database )

  ! This subroutine destroys a vector template database

    ! Dummy argument
    type (VectorTemplate_T), dimension(:), pointer :: database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: l2gpIndex, S, Status

    if ( associated(database) ) then
      do l2gpIndex = 1, SIZE(database)
         call DestroyVectorTemplateInfo ( database(l2gpIndex) )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, moduleName, 'database', s, address=addr )
    end if
  end subroutine DestroyVectorTemplateDatabase

  ! ----------------------------------  DestroyVectorTemplateInfo  -----
  subroutine DestroyVectorTemplateInfo ( VectorTemplate )

  ! This subroutine destroys a vector template created above

    ! Dummy arguments
    type (VectorTemplate_T), intent(inout) :: vectorTemplate

    ! Executable code

    call deallocate_test ( vectorTemplate%quantities, &
      & "vectorTemplate%quantities", ModuleName )

    vectorTemplate%noQuantities = 0
    vectorTemplate%totalInstances = 0
    vectorTemplate%totalElements = 0
    vectorTemplate%name = 0
    vectorTemplate%where = where_t(0,0)
  end subroutine DestroyVectorTemplateInfo

  ! -----------------------------------------  DestroyVectorValue  -----
  subroutine DestroyVectorValue ( Vector )
  ! Destroy the "values" field in all of the quantities in a vector.  This
  ! is useful when forming normal equations little-by-little.
    type(vector_T), intent(inout) :: Vector

    integer :: I

    if ( .not. associated(vector%quantities) ) return
    do i = 1, size(vector%quantities)
      call destroyVectorQuantityValue ( vector%quantities(i) )
    end do
  end subroutine DestroyVectorValue

  ! ---------------------------------------  DiffVectorQuantities  -----
  subroutine DiffVectorQuantities ( Qty1, Qty2, Name, Options )

    type (VectorValue_T), intent(in) :: QTY1, QTY2
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: OPTIONS

    if ( present(name) ) then
      if ( len_trim(name) > 0 ) then
        call output ( name )
        call output ( ', ' )
      end if
    end if
    call output ( ' Qty_Template_Name = ' )
    if ( qty1%template%name /= 0 ) then
      call display_string ( qty1%template%name )
    else
      call output ( '<none given>' )
    end if
    if ( qty1%label /= 0 ) then
      call output ( ', label = ' )
      call display_string ( qty1%label )
    else
      call output ( ' unlabeled ', advance='yes' )
    end if
    call diff ( qty1%values, 'quantity(1) values', &
      &         qty2%values, 'quantity(2) values', options=options )
  end subroutine DiffVectorQuantities

  ! ----------------------------------------------  DivideVectors  -----
  subroutine DivideVectors ( A, X, Y )
  ! Y = A / X if Y is present, else A = A / X, element-by-element.
  ! If the mask field of X's vector value is associated, only do the
  ! divide where the m_linalg bit is zero.

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: A
    type(Vector_T), intent(in) :: X
    type(Vector_T), intent(inout), optional, target :: Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    type(Vector_T), pointer :: Z
    ! Executable statements:
    if ( x%template%name /= a%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Numerator vector has different template from denominator' )
    z => a
    if ( present(y) ) then
      if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Quotient vector has different template from original' )
      z => y
    end if
    do i = 1, size(x%quantities)
      if ( associated(x%quantities(i)%mask) ) then
        where ( iand(ichar(x%quantities(i)%mask),m_linalg) == 0 )
          z%quantities(i)%values = a%quantities(i)%values / x%quantities(i)%values
        elsewhere
          z%quantities(i)%values = 0
        end where
      else
        z%quantities(i)%values = a%quantities(i)%values / x%quantities(i)%values
      end if
    end do
  end subroutine DivideVectors

  ! -------------------------------------------------  DotVectors  -----
  real(rv) function DotVectors ( X, Y ) result (Z)
  ! Compute the inner product of two vectors.

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Executable statements:
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot .DOT. vectors having different templates" )
    z = sum( x%quantities .dot. y%quantities )
  end function DotVectors

  ! -------------------------------------------  DotVectorsMasked  -----
  real(rv) function DotVectorsMasked ( X, Y ) result (Z)
  ! Compute the inner product of two vectors.  Ignore elements masked
  ! by m_linAlg in either X or Y.

    ! Dummy arguments:
    type(Vector_T), intent(in) :: X, Y
    ! Executable statements:
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot .MDOT. vectors having different templates" )
    z = sum( x%quantities .mdot. y%quantities )
  end function DotVectorsMasked

  ! --------------------------------------  DotVectorsMaybeMasked  -----
  real(rv) function DotVectorsMaybeMasked ( X, Y, UseMask ) result (Z)
  ! Compute the inner product of two vectors.  If UseMask is present and
  ! true, the masks are respected.
    type(vector_t), intent(in) :: X, Y
    logical, intent(in), optional :: UseMask
    logical MyMask
    myMask = .false.
    if ( present(useMask) ) myMask = useMask
    if ( myMask ) then
      z = x .mdot. y
    else
      z = x .dot. y
    end if
  end function DotVectorsMaybeMasked

  ! ----------------------------------------  DotVectorQuantities  -----
  elemental real(rv) function DotVectorQuantities ( X, Y ) result ( Z )
  ! Compute the inner product of two vector quantities
    type(VectorValue_T), intent(in) :: X, Y
    z = sum( x%values * y%values )
  end function DotVectorQuantities

  ! ----------------------------------  DotVectorQuantitiesMasked  -----
  elemental real(rv) function DotVectorQuantitiesMasked ( X, Y ) result ( Z )
  ! Compute the inner product of two vector quantities, respecting their
  ! masks, if any
    type(VectorValue_T), intent(in) :: X, Y
    if ( associated(x%mask) ) then
      if ( associated(y%mask) ) then
        z = sum( x%values * y%values, &
               & iand(iand(ichar(x%mask), ichar(y%mask)), m_linAlg) == 0 )
      else
        z = sum( x%values * y%values, iand(ichar(x%mask), m_linAlg) == 0 )
      end if
    else if ( associated(y%mask) ) then
      z = sum( x%values * y%values, iand(ichar(y%mask), m_linAlg) == 0 )
    else
      z = sum( x%values * y%values )
    end if
  end function DotVectorQuantitiesMasked

  ! -----------------------------  DotVectorQuantitiesMaybeMasked  -----
  elemental real(rv) function DotVectorQuantitiesMaybeMasked ( X, Y, UseMask ) &
    & result ( Z )
  ! Compute the inner product of two vector quantities, respecting their
  ! masks, if any, if UseMask is present and true
    type(VectorValue_T), intent(in) :: X, Y
    logical, intent(in), optional :: UseMask
    logical :: MyMask
    myMask = .false.
    if ( present(useMask) ) myMask = useMask
    if ( myMask ) then
      z = x .mdot. y
    else
      z = x .dot. y
    end if
  end function DotVectorQuantitiesMaybeMasked

  ! ----------------------------------------  DumpNiceMaskSummary  -----
  ! This routine tries to produce a useful human readable dump of a vector
  ! mask.
  subroutine DumpNiceMaskSummary ( qty, prefix, masksToDump )
    type(VectorValue_T), intent(in) :: qty
    character(len=*), intent(in), optional :: prefix
    integer, dimension(:), intent(in) :: masksToDump
    ! Local variables
    integer :: I, C, M                  ! Loop counters
    integer :: S0, S1                   ! Delimiters
    character, dimension ( qty%template%instanceLen ) :: MERGEDMASK ! Merge of all the masks we have
    character, dimension ( qty%template%noSurfs ) :: EXTRACTEDMASK ! Subset of merged mask for a channel
    logical, dimension ( qty%template%noSurfs ) :: THISMASK ! This particular mask flags
    logical :: VARIES                   ! Masking varies from instance to instance
    logical :: CONTINUOUS               ! Masking is in one continuous block
    character(len=2) :: CHANOFFSET      ! My prefix

    ! Executable code
    if ( .not. associated ( qty%mask ) ) then
      if ( present ( prefix ) ) call output ( prefix )
      call output ( 'This quantity has no subsetting in force.', advance='yes' )
    else
      mergedMask = qty%mask(:,1)
      varies = .false.
      do i = 2, qty%template%noInstances
        if ( any ( qty%mask(:,i) /= mergedMask ) ) then
          varies = .true.
          mergedMask = char ( iand ( ichar ( qty%mask(:,i) ), ichar ( mergedMask ) ) )
        end if
      end do
      if ( present ( prefix ) ) call output ( prefix )
      if ( varies ) then
        call output ( 'Subsetting varies from instance to instance, printing widest', advance='yes' )
      else
        call output ( 'Subsetting is constant for all instances', advance='yes' )
      end if
      do c = 1, qty%template%noChans
        chanOffset = ''
        if ( qty%template%noChans /= 1 ) then
          if ( present ( prefix ) ) call output ( prefix )
          call output ( 'Channel: ' )
          if ( size ( masksToDump ) == 1 ) then
            call output ( c )
            call output ( ' - ' )
          else
            call output ( c, advance='yes' )
            chanOffset = '  '
          end if
        end if
        extractedMask = (/ (mergedMask(c+(i-1)*qty%template%noChans), i=1, qty%template%noSurfs) /)
        ! Now have the flags for this channel, lets describe them in a helpful manner
        do m = 1, size ( masksToDump )
          thisMask = iand ( ichar ( extractedMask ), masksToDump(m) ) == 1
          call output ( chanOffset )
          ! Do a nice name for the mask
          select case ( masksToDump ( m ) )
          case ( m_ignore )
            call output ( 'Ignored: ' )
          case ( m_linAlg )
            call output ( 'Retrieved / used: ' )
          case ( m_tikhonov )
            call output ( 'Smoothed: ' )
          case ( m_fill )
            call output ( 'Fill: ' )
          case ( m_fullDerivatives )
            call output ( 'Full derivatives: ' )
          case ( m_spare )
            call output ( 'Spare mask: ' )
          end select
              
          s0 = FindFirst ( .not. thisMask )
          if ( s0 /= 0 ) then
            s1 = FindFirst ( thisMask(s0:) )
            if ( s1 == 0 ) then
              s1 = qty%template%noSurfs
            else
              s1 = s1 + s0 - 2
            end if
            continuous = count ( .not. thisMask ) == s1 - s0 + 1
            if ( s0 /= 1 .and. s1 /= qty%template%noSurfs ) then
              call outputNiceSurface ( qty%template%surfs(s0,1), qty%template%verticalCoordinate )
              call output ( ' -> ' )
              call outputNiceSurface ( qty%template%surfs(s1,1), qty%template%verticalCoordinate )
            else
              call output ( 'Everywhere' )
            end if
            if ( continuous ) then
              call output ( ' uninterrupted.', advance='yes' )
            else
              call output ( ' with gaps', advance='yes' )
            end if
          else
            call output ( ' Nowhere', advance='yes' )
          end if
        end do

      end do
    end if
  contains
    ! - - - - - - - - - - -  OutputNiceSurface
    subroutine OutputNiceSurface ( value, coordinate )
      use Intrinsic, only: L_Angle, L_Geodaltitude, L_Gph, L_None, &
        & L_Pressure, L_Theta, L_Zeta
      real(r8), intent(in) :: VALUE
      integer, intent(in) :: COORDINATE
      ! Executable code
      real(r8) :: P
      select case ( coordinate )
      case ( l_angle )
        call output ( value, format='(f0.1)' )
        call output ( ' degrees' )
      case ( l_geodAltitude, l_gph )
        call output ( value/1e3, format='(f0.2)' )
        call output ( ' km' )
      case ( l_none )
        call output ( 'No surface' )
      case ( l_pressure, l_zeta )
        if ( coordinate == l_pressure ) then
          p = value
        else
          p = 10.0**(-value)
        end if
        if ( p > 10.0 ) then
          call output ( nint ( p ) )
        else if ( p >= 1.0 ) then
          call output ( p, format='(F0.1)' )
        else if ( p >= 0.1 ) then
          call output ( p, format='(F0.2)' )
        else if ( p >= 0.01 ) then
          call output ( p, format='(F0.3)' )
        else
          call output ( p, format='(F0.5)' )
        end if
        call output ( ' hPa' )
      case ( l_theta )
        call output ( value, format='(i0)' )
        call output ( ' K' )
      end select
    end subroutine OutputNiceSurface

  end subroutine DumpNiceMaskSummary

  ! -----------------------------------------  DumpQuantityMask  -----
  ! Dump a quantity's mask with details and options set by those args:
  ! Details             effect is to print ..
  ! -------             ---------------------
  !   -1        how many locations have their bitnumth bit set 
  !   0         how many unique bitnumbers have been set and their meanings
  !   1         a table of mask bits by location
  
  ! options contain   effect
  ! ---------------   ------
  !    b[n]         set bitnum to n (default was 0)
  !     l           collapse locations to just instances, showing
  !                   only if bitnumth bit set at all heights, channels
  !     L           collapse locations to just instances, showing
  !                   only if bit NOT set at all heights, channels
  !
  subroutine DumpQuantityMask ( VECTORQUANTITY, DETAILS, OPTIONS )
    use MLSStringLists, only: ExpandStringRange, OptionDetail
    use MLSStrings, only: ReadNumsFromChars
    type (VectorValue_T), intent(in)       :: VectorQuantity
    integer, intent(in), optional          :: DETAILS ! if < 0 just dump summary
    character(len=*), optional, intent(in) :: options

    ! Local variables
    character(len=32) :: bitCh
    integer :: bitNum                   ! bit number; e.g., 0 is m_linAlg
    integer, dimension(18) :: bitNums
    integer :: c                        ! Channel index
    integer :: i                        ! Instance index
    integer :: j                        ! Element index
    integer :: k                        ! bitNums index
    integer :: myDetails
    character(len=32) :: myOptions
    integer :: n
    integer :: nBitNums
    integer :: nUnique
    integer :: s                        ! Surface index
    integer, dimension(10000) :: uniqueVals
    integer :: w                        ! Line width used so far
    ! Executable
    nUnique = 0
    uniqueVals = 0
    myDetails = 0
    if ( present(details) ) myDetails = details
    myOptions = ' '
    if ( present(options) ) myOptions = options
    bitNum = 0
    nBitNums = -1
    if ( index(myOptions, 'b') > 0 ) then
      bitCh = optionDetail( myOptions, single_option='b' )
      if ( index( bitCh, ',' ) > 0 .or. index( bitCh, '-' ) > 0 ) then
        call ExpandStringRange ( bitCh, bitNums, nBitNums )
      else
        call readNumsFromChars( bitCh, bitnum )
      end if
    end if
    if ( nBitNums < 0 ) then
      nBitNums = 1
      bitNums(1) = bitNum
    end if
    call output ( 'Quantity ' )
    call display_string ( vectorQuantity%template%name )

    if ( .not. vectorQuantity%template%regular ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Unable to dump mask for irregular quantities' )

    if ( .not. associated ( vectorQuantity%mask ) ) then
      call output ( ' has no mask.', advance='yes' )
      return
    else if( myDetails == 0 ) then
      call dump ( ichar(vectorQuantity%mask), name='  Mask =', &
        & format='(z3.2)', width = 20 )
      call FindUnique( &
      & reshape( ichar(vectorQuantity%mask), (/ size(vectorQuantity%mask,1)*size(vectorQuantity%mask,2) /) ), &
      & uniqueVals, nUnique )
      call outputNamedValue(' Number unique values', nUnique )
      do i=1, nUnique
        call DumpBitNames( uniqueVals(i), MaskBitNames )
      end do
      return
    else if ( index(myOptions, 'l') < 1 .and. index(myOptions, 'L') < 1 ) then
      call newLine
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
              call newLine
              call output ( '      ' )
              w = 6
            end if
            w = w + 3
            j = j + 1
            call output ( ichar(vectorQuantity%mask(j,i)), &
              & format='(z3.2)' )
          end do
          call newLine
        end do                        ! Surface loop
      end do                          ! Instance loop
      return
    end if
    ! These are the cases where we show results broken down by bit numbers
    do k=1, nBitNums
      nUnique = 0
      uniqueVals = 0
      bitnum = bitNums(k)
      if ( nBitNums > 1 ) then
        call newLine
        call outputnamedvalue( 'bit number', bitnum )
        call DumpBitNames( 2**bitnum, MaskBitNames, oneLine=.true. )
      end if
      if( myDetails < 0 ) then
        n = count( isBitSet( ichar(vectorQuantity%mask), bitNum ) )
        call outputNamedValue( 'Num mask bits set: ', n )
      else if ( index(myOptions, 'l') > 0 ) then
        ! Collapse locations, showing where entirely masked channels and surfaces
        do i = 1, vectorQuantity%template%noInstances
          if ( any( &
            & .not. isBitSet( ichar(vectorQuantity%mask(:,i)), bitNum ) &
            & ) ) cycle
          nUnique = nUnique + 1
          uniqueVals ( nUnique ) = i
        end do
        if ( nUnique > 0 ) then
          call dump( uniqueVals(1:nUnique), ': Locations entirely masked' )
        else
          call output( ': No locations entirely masked', advance = 'yes' )
        end if
      else if ( index(myOptions, 'L') > 0 ) then
        ! Collapse locations, showing where NOT masked
        do i = 1, vectorQuantity%template%noInstances
          if ( all( &
            & isBitSet( ichar(vectorQuantity%mask(:,i)), bitNum ) &
            & ) ) cycle
          nUnique = nUnique + 1
          uniqueVals ( nUnique ) = i
        end do
        if ( nUnique > 0 ) then
          call dump( uniqueVals(1:nUnique), ': Locations not entirely masked' )
        else
          call output( ': All locations entirely masked', advance = 'yes' )
        end if
      end if
    end do
  end subroutine DumpQuantityMask

  ! ---------------------------------------------  DumpVectorMask  -----
  subroutine DumpVectorMask ( VECTOR, DETAILS, OPTIONS )
    type (Vector_T), intent(in) :: VECTOR
    integer, intent(in), optional :: DETAILS ! if < 0 just dump summary
    character(len=*), optional, intent(in) :: options

    ! Local variables
    integer :: q                        ! Quantity index

    ! Executable code
    call output ( 'Dumping mask for vector ' )
    call display_string ( vector%name, advance='yes' )

    do q = 1, size(vector%quantities)
      call dumpMask ( vector%quantities(q), details, options )
    end do                              ! Loop over quantities
  end subroutine DumpVectorMask
  ! --------------------------------------------  DumpVectorNorms  -----
  subroutine DumpVectorNorms ( Vector, Level, Name, UseMask )
    type (Vector_T), intent(in) :: Vector
    integer, intent(in) :: Level ! <=0 => Whole vector norm, 
                                 ! 1 => quantity norms,
                                 ! >= 2 => quantity norms with quantity names
    character(*), intent(in), optional :: Name
    logical, intent(in), optional :: UseMask
    integer :: Q
    select case ( level )
    case ( :0 )
      if ( present(name) ) call output ( name )
      call output ( l2norm ( vector, useMask ) )
    case ( 1 )
      call dump ( l2norm ( vector%quantities, useMask ), name )
    case ( 2 )
      if ( present(name) ) call output ( name, advance='yes' )
      do q = 1, size(vector%quantities)
        call output ( q, format='(i4)', after='# ' )
        call display_string ( vector%quantities(q)%template%name )
        call output ( l2norm(vector%quantities(q),useMask), before=': ', &
          & advance='yes' )
      end do
    end select
  end subroutine DumpVectorNorms

  ! ------------------------------------------------  Dump_Vector  -----
  subroutine Dump_Vector ( VECTOR, DETAILS, OPTIONS, NAME, &
    & QUANTITYTYPES, INSTRUMENTMODULES, SIGNAL_IDS, &
    & COHERENT, STACKED, REGULAR, MINORFRAME, MAJORFRAME, &
    & THENDITCHAFTERDUMP, CLEAN )

    use Lexer_Core, only: Print_Source

    ! dump quantities in vector according to whether they match
    ! all of the optional args: name, ..,majorframe
    ! if thenditchafterdump is present and TRUE,
    ! dump only first matching quantity
    type(Vector_T), intent(in) :: VECTOR
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump quantity values
    !                                        ! -1 Skip quantity details beyond names
    !                                        ! -2 Skip all quantity details but
                                             !    print a size summary
                                             ! <= -3 => no output
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: OPTIONS
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
    logical, intent(in), optional                :: CLEAN

    ! Local parameters
    logical :: dumpThisQty
    integer :: J    ! Loop inductor, subscript
    integer :: MyDetails
    logical :: myditchafterdump
    character(len=8) :: myOptions
    integer :: TotalSize
    
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( myDetails <= -3 ) return

    if ( present(thenditchafterdump) ) then
      myditchafterdump = thenditchafterdump
    else
      myditchafterdump = .false.
    end if

    totalSize = 0

    if ( present(name) ) then
      call output ( name ); call output ( ', ' )
    end if
    if ( vector%name /= 0 ) then
      call output ( 'Name = ' )
      call display_string ( vector%name )
    end if
    if ( vector%where%source /= 0 ) then
      call output ( ' created at ' )
      call print_source ( vector%where )
    end if
    if ( vector%template%name /= 0 ) then
      call output ( ' Template_Name = ' )
      call display_string ( vector%template%name )
    end if
    call newline
    myoptions = ' '
    if ( present(clean) ) then
      if ( clean ) myoptions = 'c'
    end if
    if ( present(options) ) myOptions = trim(myOptions) // options
    do j = 1, size(vector%quantities)
      dumpThisQty = myDetails > -2
      if ( associated(vector%quantities(j)%values) ) &
        & totalSize = totalSize + size(vector%quantities(j)%values)
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
        call output ( j, 4, after="~" )
        call dump ( vector%quantities(j), details, options=myOptions )
        if ( myditchafterdump ) return
      end if
    end do ! j
    if ( myDetails == -2 ) then
      call output ( size(vector%quantities), before='      having ' )
      call output ( totalSize, before=' quantities and ' )
      call output ( ' elements.', advance='yes' )
    end if
  end subroutine Dump_Vector

  ! -----------------------------------------------  Dump_Vectors  -----
  subroutine Dump_Vectors ( VECTORS, DETAILS, OPTIONS, NAME, &
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
    !                                        ! -3 Just summarize the database
    !                                        ! >0 Do dump quantity values
    !                                        ! Default 1
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: OPTIONS
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
    integer :: MyDetails
    logical :: MyDitchAfterDump
    integer :: TotalSize ! of vector quantities

    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( present(thenditchafterdump) ) then
      myditchafterdump = thenditchafterdump
    else
      myditchafterdump = .false.
    end if
    if ( size(vectors) > 1 ) &
      & call output ( size(vectors), before='VECTORS: SIZE = ', advance='yes' )

    totalSize = 0

    do i = 1, size(vectors)
      ! Presume do not need to dump vector; hence preset to FALSE -- 
      ! becomes TRUE if wish to dump one or more quantities
      dumpThisVector = .false.
      if ( .not. associated(vectors(i)%quantities) .and. myDetails > -3 ) then
        call output ( '(entry  ', advance='no' )
        call output ( i, advance='no' )
        call output ( '  in the vector database had been destroyed)  ', &
        & advance='yes' )
        cycle
      end if
      do j=1, size(vectors(i)%quantities)
        if ( associated(vectors(i)%quantities(j)%values) ) &
          totalSize = totalSize + size(vectors(i)%quantities(j)%values)
        ! Presume need to dump quantity; hence preset to TRUE --
        ! becomes FALSE if fails to match a requirement
        dumpThisQty = myDetails > -3
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
      end do
      if ( dumpThisVector ) then
        call output ( i, 4 )
        call output ( ': ' )
        call dump_vector ( vectors(i), details, options, name, &
        & quantitytypes, instrumentmodules, signal_ids, &
        & coherent, stacked, regular, minorframe, majorframe, &
        & thenditchafterdump )
        if ( myditchafterdump ) return
      end if
    end do ! i
    call output ( size(vectors), before='Vectors database has ' )
    call output ( totalSize, before=' vectors with ' )
    call output ( ' elements.', advance='yes' )
  end subroutine Dump_Vectors

  ! ---------------------------------------  Dump_Vector_Quantity  -----
  subroutine Dump_Vector_Quantity ( Qty, Details, Name, Vector, Options )

    use Pointer_Rank_Remapping, only: Remap

    type (VectorValue_T), intent(in) :: Qty
    integer, intent(in), optional :: Details ! <0  => Name only
    !                                        ! <-1 => No newline after name
    !                                        ! =0  => No values
    !                                        ! >0  => Dump quantity values
    !                                        ! >1  => Dump template with details-1
    !                                        ! Default 1
    character(len=*), intent(in), optional :: Name
    type (Vector_T), intent(in), optional  :: Vector ! Only to get its name
    character(len=*), intent(in), optional :: Options ! E.g., '-sb'
    ! If options is present and
    ! contains        dump        but skip
    !    1            values     template data
    !    2            values     most template data
    !    3        values and template data
    ! Internal variables
    logical, parameter :: deeBug = .false.  ! Should this also affect Details?
    logical :: Dot ! Use vector.quantity notation
    integer :: i
    character(len=32) :: oldInfo
    integer :: myDetails
    character(len=8) :: myOptions
    integer :: nUnique
    integer, dimension(1000) :: uniqueVals
    real(rv), pointer :: value4(:,:,:,:) ! Remapping of Qty%values

    ! Executable
    myDetails = 1
    if ( present(details) ) myDetails = details
    myOptions = ' '
    if ( present(options) ) myOptions = options
    oldInfo = MLSMessageConfig%Info
    if ( present(name) ) then
      MLSMessageConfig%Info = name
      call output ( name ); call output ( ', ' )
    else if ( qty%template%name /= 0 ) then
      call get_string ( qty%template%name, MLSMessageConfig%Info )
    end if
    if ( .not. index(myOptions, '1') > 0 ) then
      dot = .false.
      if ( present(vector) ) dot = vector%name /= 0
      if ( dot ) then
        call output ( ' Vector quantity name = ' )
        if ( vector%name /= 0 ) then
          call display_string ( vector%name )
        else
          call output ( '<none given>' )
        end if
        if ( qty%template%name /= 0 ) then
          call display_string ( qty%template%name, before='.' )
        else
          call output ( '.<none given>' )
        end if
      else
        call output ( ' Qty_Template_Name = ' )
        if ( qty%template%name /= 0 ) then
          call display_string ( qty%template%name )
        else
          call output ( '<none given>' )
        end if
      end if
      if ( qty%label /= 0 ) then
        call output ( ', label = ' )
        call display_string ( qty%label )
      else
        call output ( ' unlabeled ', advance='yes' )
      end if
      if ( myDetails < -1 ) then
        MLSMessageConfig%Info = oldInfo
        return
      end if
      call newLine
      if ( myDetails < 0 ) then
        MLSMessageConfig%Info = oldInfo
        return
      end if
    else if ( qty%template%name /= 0 ) then
      call display_string ( qty%template%name, before='.' )
    end if
    if ( .not. index(myOptions, '1') > 0 &
      & .and. .not. index(myOptions, '2') > 0 ) then
      if ( qty%template%quantityType == l_vmr ) then
        call output ( '    molecule: ')
        if ( qty%template%molecule < 1 ) then
          call output ( '    (no database entry for this quantity) ' )
        else
          call display_string ( lit_indices(qty%template%molecule) )
        end if
      else if ( qty%template%quantityType > 0 ) then
        call output( '    Quantity type: ', advance='no' )
        call display_string ( lit_indices(qty%template%quantityType), advance='yes' )
      end if
      call output ( qty%template%noChans, before='    noChans = ' )
      call output ( qty%template%noSurfs, before=' noSurfs = ' )
      call output ( qty%template%noInstances, before=' noInstances = ')
      call output ( qty%template%noCrossTrack, before=' noCrossTrack = ' )
      call output ( qty%template%instanceLen, before=' instanceLen = ', advance='yes')
      call output ( '    signal: ')
      if ( qty%template%signal < 1 ) then
        call output ( '    (no database entry for this quantity) ', advance='yes')
      else if ( signals(qty%template%signal)%name < 1 ) then
        call output ( '    (no name in the database for this quantity) ', advance='yes')
      else
        call display_string ( signals(qty%template%signal)%name, advance='yes' )
      end if
      call output ( '    instrumentmodule: ')
      if ( qty%template%instrumentModule < 1 ) then
        call output ( '    (no database entry for this quantity)', advance='yes' )
      else
        call display_string ( modules(qty%template%instrumentModule)%name, advance='yes' )
      end if
      call output ( '    ' )
      if ( .not. qty%template%coherent ) call output ( 'in' )
      call output ( 'coherent ' )
      if ( .not. qty%template%stacked ) call output ( 'non' )
      call output ( 'stacked ' )
      if ( .not. qty%template%regular ) call output ( 'ir' )
      call output ( 'regular ' )
      if ( qty%template%logBasis ) then
        call output ('log-')
      else
        call output ('linear-')
      end if
      call output ('basis ' )
      call output ( trim(merge('   ','non',qty%template%minorFrame)) // &
        & 'minorFrame' )
      call output ( ' ' // trim(merge('   ','non',qty%template%majorFrame)) // &
      & 'majorFrame', advance='yes' )
      if ( size(qty%values) > 0 ) then
        call output ( '    values array size is ' )
        if ( qty%template%noChans > 1 ) then
          call output ( qty%template%noChans )
          call output ( 'x' )
        end if
        call output ( qty%template%noSurfs )
        call output ( qty%template%noInstances, before='x' )
        if ( qty%template%noCrossTrack > 1 ) &
          & call output ( qty%template%noCrossTrack, before='x' )
        call output ( size(qty%values), before=' = ' )
      else
        call output ( '    values array size is 0' )
      end if
    end if
    if ( .not. associated(qty%values) ) then
      call output( 'values array is not associated', advance='yes' )
    else if ( size(qty%values) < 1 ) then
      call output( 'values is a 0-size array', advance='yes' )
    else if ( myDetails > 0 ) then
      call newLine
      call remap ( qty%values, value4, &
        & [ qty%template%noChans, qty%template%noSurfs, &
        &   qty%template%noInstances, qty%template%noCrossTrack ] )
    ! value4(1:qty%template%noChans,1:qty%template%noSurfs, &
    !       &1:qty%template%noInstances,1:qty%template%noCrossTrack) => qty%values
      call dump ( value4, '  Elements = ', options=options )
      if ( associated(qty%mask) ) then
        call dump ( ichar(qty%mask), name='  Mask(hex) =', &
          & format='(z3.2)', width = 20 )
!           call dumpQuantityMask ( qty )
        call FindUnique( &
        & reshape( ichar(qty%mask), (/ size(qty%mask,1)*size(qty%mask,2) /) ), &
        & uniqueVals, nUnique )
        call outputNamedValue(' Number unique values', nUnique )
        do i=1, nUnique
          call DumpBitNames( uniqueVals(i), MaskBitNames )
        end do
      else
        call output ( '      Without mask', advance='yes' )
      end if
    else
      call output ( ', with' )
      if ( .not. associated(qty%values) ) call output ( 'out' )
      call output ( ' values, with' )
      if ( .not. associated(qty%mask ) ) call output ( 'out' )
      call output ( ' mask', advance='yes' )
    end if
    if ( associated(qty%BinNumber) ) then
      call Dump( qty%BinNumber, 'Bin Numbers' )
    elseif ( deebug ) then
      call output( 'Qty BinNumber not associated', advance='yes' )
    endif
    if ( associated(qty%MAF) ) then
      call Dump( qty%MAF, 'MAFs' )
    elseif ( deebug ) then
      call output( 'Qty MAF not associated', advance='yes' )
    endif
    if ( myDetails > 1 ) call dump ( qty%template, details=myDetails-1 )
    MLSMessageConfig%Info = oldInfo
  end subroutine Dump_Vector_Quantity

  ! ---------------------------------------  Dump_Vector_Template  -----
  subroutine Dump_Vector_Template ( VECTOR_TEMPLATE, DETAILS, QUANTITIES )
    use QuantityTemplates, only: Dump
    type(VectorTemplate_T), intent(in) :: VECTOR_TEMPLATE
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
                                             ! > 0  => Do dump arrays
                                             ! Default 1
    type(quantityTemplate_T), intent(in), optional :: QUANTITIES(:)
                                             ! If DETAILS > 0 and QUANTITIES
                                             ! is present, dump them with level
                                             ! details-1.
    integer :: I, MyDetails
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( vector_template%name /= 0 ) then
      call output ( ' Name = ' )
      call display_string ( vector_template%name )
    end if
    call output ( ' NoQuantities = ' )
    call output ( vector_template%noQuantities )
    call output ( ' TotalInstances = ' )
    call output ( vector_template%totalInstances )
    call output ( ' TotalElements = ' )
    call output ( vector_template%totalElements, advance='yes' )
    if ( myDetails > 0 ) then
      call dump ( vector_template%quantities, '      Quantities = ' )
      if ( present(quantities) ) then
        do i = 1, size(vector_template%quantities)
          call output ( vector_template%quantities(i), 4 )
          call output ( ':' )
          call dump ( quantities(vector_template%quantities(i)), details-1 )
        end do
      end if
    else if ( .not. associated(vector_template%quantities) ) then
       call output ( 'vector quantities not associated', advance='yes' )
    else
       call outputNamedValue( 'size(vector_template%quantities)', size(vector_template%quantities) )
    end if
  end subroutine Dump_Vector_Template

  ! --------------------------------------  Dump_Vector_Templates  -----
  subroutine Dump_Vector_Templates ( VECTOR_TEMPLATES, DETAILS )
    type(VectorTemplate_T), intent(in) :: VECTOR_TEMPLATES(:)
    integer, intent(in), optional :: DETAILS ! <= 0 => Don't dump arrays
    !                                        ! > 0  => Do dump arrays
    !                                        ! Default 1
    integer :: I
    call output ( 'VECTOR_TEMPLATES: SIZE = ' )
    call output ( size(vector_templates), advance='yes' )
    do i = 1, size(vector_templates)
      call output ( i, 4 )
      call output ( ':' )
      call dump_vector_template ( vector_templates(i), details=details )
    end do
  end subroutine Dump_Vector_Templates

  ! ------------------------------------------  GatherVectorQuantity  -----
  function GatherVectorQuantity ( quantity, start, count, stride, block )

  ! This function returns a pointer to the information about one quantity
  ! within a vector.

    use Pointer_Rank_Remapping, only: Remap

    ! Dummy arguments
    type(VectorValue_T), pointer :: Quantity
    integer, dimension(:), intent(in) :: start
    integer, dimension(:), intent(in) :: count
    integer, dimension(:), intent(in) :: stride
    integer, dimension(:), intent(in) :: block

    ! Result
    type(VectorValue_T) :: GatherVectorQuantity

    ! Internal
    real(rv), pointer :: GV1(:), QV1(:), GV3(:,:,:), QV3(:,:,:)

    ! Executable
    call CloneVectorQuantity( GatherVectorQuantity, quantity )
    call DestroyVectorQuantityValue( GatherVectorQuantity, destroyMask=.true., &
      & destroyTemplate=.false. )
    select case( size(count) )
    case (1)
    ! qv1(1:size(quantity%values)) => quantity%values
      call remap ( quantity%values, qv1, size(quantity%values) )
      call ExtractArray ( gv1, qv1, start, count, stride, block, options='-a' )
      GatherVectorQuantity%values(1:quantity%template%noChans*quantity%template%noSurfs, &
        & 1:quantity%template%noInstances*quantity%template%noCrossTrack) => gv1
    case (2)
      call ExtractArray ( GatherVectorQuantity%values, quantity%values, &
        & start, count, stride, block, options='-a' )
    case (3)
    ! qv3(1:quantity%template%noChans,1:quantity%template%noSurfs, &
    !    &1:quantity%template%noInstances*quantity%template%noCrossTrack) => &
    !    & quantity%values
      call remap ( quantity%values, qv3, [ quantity%template%noChans, &
        & quantity%template%noSurfs, &
        & quantity%template%noInstances*quantity%template%noCrossTrack ] )
      call ExtractArray ( gv3, qv3, start, count, stride, block, options='-a' )
    ! GatherVectorQuantity%values(1:quantity%template%noChans*quantity%template%noSurfs, &
    !  1:quantity%template%noInstances*quantity%template%noCrossTrack) => gv3
      call remap ( gv3, GatherVectorQuantity%values, &
        & [ quantity%template%noChans*quantity%template%noSurfs, &
        &   quantity%template%noInstances*quantity%template%noCrossTrack ] )
    case default
      call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "GatherVectorQuantity can handle only count sized 1, 2, or 3" )
    end select
  end function GatherVectorQuantity

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
    & molecule, instrumentModule, supportedInstrumentModule, radiometer, reflector, signal, &
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
    integer, intent(in),  optional :: SUPPORTEDINSTRUMENTMODULE ! Another instrument module index
    integer, intent(in),  optional :: RADIOMETER   ! Radiometer index
    integer, intent(in),  optional :: REFLECTOR   ! Reflector literal
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
    if ( present(noError)) myNoError = noError

    GetVectorQuantityByType => NULL()

    ! Look in the first vector
    index = GetVectorQuantityIndexByType ( vector, &
      & quantityType, molecule, instrumentModule, supportedInstrumentModule, radiometer, reflector, signal, &
      &   sideband, noError = present(otherVector) .or. myNoError)
    if ( index /= 0 ) then
      if ( present (foundInFirst) ) foundInFirst = .true.
      GetVectorQuantityByType => vector%quantities(index)
    else
      ! Can only get here if not found in first vector and noError or other
      ! vector
      if ( present (otherVector) ) then
        index = GetVectorQuantityIndexByType ( otherVector, &
          &  quantityType, molecule, instrumentModule, supportedInstrumentModule, radiometer, reflector, signal, &
          &  sideband, noError=myNoError )
        if ( present (foundInFirst) ) foundInFirst = .false.
        if ( index /= 0 ) &
          & GetVectorQuantityByType => otherVector%quantities( index )
      end if
    end if
  end function GetVectorQuantityByType

  ! --------------------------------  GetVectorQtyByTemplateIndex  -----
  ! Given a vector and an index into the quantity templates, find quantity
  ! with matching template within vector.
  function GetVectorQtyByTemplateIndex ( vector, quantityIndex, indexInVector )
    ! The vector where the quantity resides
    type (vector_T), intent(in) :: vector
    ! The index of the quantity in the quantity database; this index
    ! is returned by AddQuantityTemplateToDatabase
    integer, intent(in) :: quantityIndex
    integer, intent(out), optional :: indexInVector
    ! Result
    type (VectorValue_T), pointer :: GetVectorQtyByTemplateIndex

    ! Local variables
    integer :: myIndexInVector
    character(len=132) :: MSG           ! An error message

    ! Executable code
    GetVectorQtyByTemplateIndex => NULL()
    if ( .not. associated ( vector%quantities ) ) then
      msg = "Reference to the vector "
      call get_string ( vector%name, msg(len_trim(msg)+2:len(msg)), strip=.true. )
      msg = trim ( msg ) // " which has been destroyed" 
      call MLSMessage ( MLSMSG_Error, ModuleName, trim(msg) )
    end if
    myIndexInVector = FindFirst ( vector%template%quantities, quantityIndex )
    if ( myIndexInVector > 0 .and. &
      & myIndexInVector <= ubound(vector%quantities, 1) ) &
      & GetVectorQtyByTemplateIndex => &
      &   vector%quantities(myIndexInVector)
    if ( myIndexInVector > size(vector%quantities) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'index returned by FindFirst too big' )
    if ( present ( indexInVector ) ) indexInVector = myIndexInVector
  end function GetVectorQtyByTemplateIndex

  ! -------------------------------  GetVectorQuantityIndexByName  -----
  ! This family of subroutines returns a quantity's index within a vector
  ! by its quantity name
  function GetVectorQuantityIndexByName_sbr ( VECTOR, QUANTITYNAME, &
    & NOERR ) result( index )

  ! Given a quantity name's sub-rosa index, this function returns the index
  ! of the quantity within the vector matching that name.

    ! Dummy arguments
    type (Vector_T), intent(in)   :: VECTOR
    integer, intent(in)           :: QUANTITYNAME ! Quantity name sub-rosa index.
    logical, intent(in), optional :: NOERR        ! No error if present and true.
                                                  ! Return -1 if no such quantity.
    integer                       :: INDEX
    ! Local variables
    character(len=127) :: MSG
    integer :: Search

    ! Executable code
    do search = 1, size(vector%quantities)
      if ( quantityName == vector%quantities(search)%template%name ) then
        index = search
        return
      end if
    end do
    if ( present(noErr) ) then
      index = -1
      if ( noErr ) return
    end if
    call get_string ( quantityName, msg )
    msg(string_length(quantityName)+2:) = 'is not a quantity in vector'
    call get_string ( vector%name, msg(len_trim(msg)+2:) )
    call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )

  end function GetVectorQuantityIndexByName_sbr

  function GetVectorQuantityIndexByName_char ( VECTOR, QUANTITYNAME, &
    & NOERR ) result( index )
    use MLSStrings, only: Lowercase

  ! Given a quantity name, this function returns the index
  ! of the quantity within the vector matching that name.

    ! Dummy arguments
    type (Vector_T), intent(in)   :: VECTOR
    character(len=*), intent(in)  :: QUANTITYNAME ! Quantity name 
    logical, intent(in), optional :: NOERR        ! No error if present and true.
                                                  ! Return -1 if no such quantity.
    integer                       :: INDEX

    ! Local variables
    character(len=127) :: MSG
    integer :: Search

    ! Executable code
    ! call output( trim(quantityName), advance='yes' )
    do search = 1, size(vector%quantities)
      call get_string ( vector%quantities(search)%template%name, msg, strip=.true. )
      ! call output( trim(msg), advance='yes' )
      if ( lowercase(quantityName) == lowercase(msg) ) then
        index = search
        return
      end if
    end do
    if ( present(noErr) ) then
      index = -1
      if ( noErr ) return
    end if
    msg = trim(quantityName) // ' is not a quantity in vector'
    call get_string ( vector%name, msg(len_trim(msg)+2:) )
    call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )

  end function GetVectorQuantityIndexByName_char

  ! -------------------------------  GetVectorQuantityIndexByType  -----
  integer function GetVectorQuantityIndexByType ( vector, quantityType, &
    & molecule, instrumentModule, supportedInstrumentModule, radiometer, &
    & reflector, signal, sideband, noError )

  ! Given a quantity type index (l_...), this function returns the index
  ! of the first quantity within the vector that has that type.  If
  ! molecule and/or radiometer are supplied, the quantity that has the
  ! specified type, as well as the specified molecule and/or radiometer
  ! index, is returned.

    use MLSSignals_m, only: GetRadiometerName
    use Molecules, only: IsExtinction

    ! Dummy arguments
    type (Vector_T), intent(in) :: VECTOR
    integer, intent(in) :: QUANTITYTYPE ! Quantity type index (l_...)
    integer, intent(in), optional :: MOLECULE     ! Molecule index (l_...)
    integer, intent(in), optional :: INSTRUMENTMODULE ! Module index
    integer, intent(in), optional :: SUPPORTEDINSTRUMENTMODULE ! Another module index
    integer, intent(in), optional :: RADIOMETER   ! Radiometer index
    integer, intent(in), optional :: REFLECTOR   ! Reflector literal
    integer, intent(in), optional :: SIGNAL       ! Signal Index
    integer, intent(in), optional :: SIDEBAND ! -1, 0, +1
    logical, intent(in), optional :: NOERROR ! Don't give error if not found
    type (QuantityTemplate_T), pointer :: QT

    ! Local variables
    character(len=127) :: MSG
    integer :: SEARCH
    logical :: MYNOERROR

    myNoError = .false.
    if ( present(noError)) myNoError = noError

    ! Executable code
    do search = 1, size(vector%quantities)
      qt => vector%quantities(search)%template
      if ( quantityType == qt%quantityType ) then
        if ( present(molecule) ) then
          if ( qt%molecule /= molecule ) cycle
        end if
        if ( present(instrumentModule) ) then
          if ( qt%instrumentModule /= instrumentModule ) cycle
        end if
        if ( present(supportedInstrumentModule) ) then
          if ( qt%instrumentModule == 0 ) cycle
          if ( modules ( qt%instrumentModule )%supportedModule /= supportedInstrumentModule ) cycle
        end if
        if ( present(radiometer) ) then
          ! Somewhat trickier than one might think this one.  We don't
          ! want to include quantities that have a radiometer purely by virtue
          ! of having a signal
          if ( qt%signal /= 0 ) cycle
          ! We can be a little lenient with vmr here.
          if ( quantityType == l_vmr ) then
            if ( .not. present ( molecule ) ) call MLSMessage ( &
              & MLSMSG_Error, ModuleName, "Requests for vmrs must have molecules" )
            if ( radiometer /= qt%radiometer .and. isExtinction(molecule) ) cycle
          else
            if ( radiometer /= qt%radiometer ) cycle
          end if
        end if
        if ( present(reflector) ) then
          if ( qt%reflector /= reflector ) cycle
        end if
        if ( present(signal) ) then
          if ( qt%signal /= signal ) cycle
        end if
        if ( present(sideband) ) then
          if ( qt%sideband /= sideband ) cycle
        end if
        GetVectorQuantityIndexByType = search
        return
       end if
    end do

    ! Not found, perhaps generate an error
    if ( myNoError ) then
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

      if ( present ( molecule ) ) then
        msg = trim(msg) // ' for molecule'
        call get_string ( lit_indices(molecule), msg(len_trim(msg)+2:))
      end if

      if ( present ( radiometer ) ) then
        msg = trim(msg) // ' for radiometer'
        call getradiometerName ( radiometer, msg(len_trim(msg)+2:))
      end if

      if ( present ( reflector ) ) then
        msg = trim(msg) // ' for reflector'
        call get_string ( lit_indices(reflector), msg(len_trim(msg)+2:))
      end if

      if ( present ( instrumentModule ) ) then
        msg = trim(msg) // ' for instrument module'
        call get_string ( modules(instrumentModule)%name, msg(len_trim(msg)+2:))
      end if

      if ( present ( signal ) ) then
        msg = trim(msg) // ' for signal'
        call GetSignalName ( signal, msg(len_trim(msg)+2:), sideband=sideband )
      end if

      call dump( vector, details=0 ) ! It's not the values causing us to crash
      call MLSMessage ( MLSMSG_Error, ModuleName, msg(:len_trim(msg)) )
    end if

  end function GetVectorQuantityIndexByType

  ! ---------------------------------------- InflateVectorDatabase -----
  integer function InflateVectorDatabase ( database, extra )
    ! Make a vector database bigger by extra
    ! Return index of first new element

    ! Dummy arguments
    type (Vector_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type (Vector_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateVectorDatabase = firstNewItem
  end function InflateVectorDatabase

  ! -------------------------------- InflateVectorTemplateDatabase -----
  integer function InflateVectorTemplateDatabase ( database, extra )
    ! Make a vector template database bigger by extra
    ! Return index of first new element

    ! Dummy arguments
    type (VectorTemplate_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type (VectorTemplate_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateVectorTemplateDatabase = firstNewItem
  end function InflateVectorTemplateDatabase

  ! ------------------------------------------  IsVectorQtyMasked  -----
  logical function IsVectorQtyMasked ( vectorQty, Row, Column, What )

  ! Is the mask for VectorQty set for address (Row, Column) ?
  ! If What is present, look at the bits of the mask specified by the union
  ! of the nonzero bits of What.  Otherwise, look at the M_LinAlg bit.
  
  ! Formal args
    type (VectorValue_T), intent(in) :: vectorQty
    integer, intent(in) ::              ROW
    integer, intent(in) ::              COLUMN
    integer, intent(in) ::              WHAT

    isVectorQtyMasked = .false.
    if ( .not. associated(vectorQty%mask)) return
    isVectorQtyMasked = iand(ichar(vectorQty%mask(row, column)), What) /= 0

  end function IsVectorQtyMasked

  ! ----------------------------------------------------  L2NormQ  -----
  elemental real(rv) function L2NormQ ( Qty, UseMask )
    ! Compute the L2Norm of a vector quantity, respecting its mask if
    ! useMask is present and true.
    type(vectorValue_t), intent(in), optional :: Qty
    logical, intent(in), optional :: UseMask
    l2NormQ = 0.0
    if ( present(qty) ) &
      & l2NormQ = sqrt( dotVectorQuantitiesMaybeMasked ( qty, qty, useMask ) )
  end function L2Normq

  ! ----------------------------------------------------  L2NormV  -----
  real(rv) function L2NormV ( Vector, UseMask )
    ! Compute the L2Norm of a vector, respecting its mask if useMask is
    ! present and true.
    type(vector_t), intent(in) :: Vector
    logical, intent(in), optional :: UseMask
    l2NormV = 0.0
    if ( associated(vector%quantities) ) &
      & l2normV = sqrt( dotVectorsMaybeMasked ( vector, vector, useMask ) )
  end function L2NormV

  ! ----------------------------------------------  MaskVectorQty  -----
  subroutine MaskVectorQty ( vectorQty, Row, Column, What )

  ! Set bits of the mask for VectorQty(Row,Column); meaning
  ! If set, don't use vectorQty%values(Row, Column)
  ! Otherwise, go ahead.  If What is present, set the bits in
  ! mask indicated by What.  Otherwise, set the M_LinAlg bit.
  
  ! Note that 
  ! (1) if mask is not yet it associated, it will be created
  ! (2) If mask already associated, any bits in it will not be reset
  !    (i.e., we will be "or"-ing bits)
  
  ! Formal args
    type (VectorValue_T), intent(inout) :: vectorQty
    integer, intent(in) ::              ROW
    integer, intent(in) ::              COLUMN
    integer, intent(in) ::              WHAT

    if ( .not. associated(vectorQty%mask)) call CreateMask ( vectorQty)
    vectorQty%mask(row, column) = &
      & char(ior( ichar(vectorQty%mask(row, column)), What ) )

  end subroutine MaskVectorQty

  ! -----------------------------------------  MoveVectorQuantity  -----
  subroutine MoveVectorQuantity ( From, To )
    ! Deallocate the Values and Mask fields in To
    ! Move the Values and Mask fields from From to To using pointer assignment.
    ! Nullify the Values and Mask fields in From
    type(vectorValue_t), intent(inout) :: From, To
    if ( from%template%name /= to%template%name ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, 'From and To vectors have different templates' )
    call destroyVectorQuantityValue ( to, destroyMask=.true. )
    to%value1 => from%value1
    to%mask1 => from%mask1
    call remapVectorValue ( to )
    call remapVectorMask ( to )
    nullify ( from%value1, from%mask1 ) ! Don't deallocate during destroy!
    nullify ( from%values, from%mask )  ! Don't deallocate during destroy!
    nullify ( from%value3, from%mask3 ) ! Don't deallocate during destroy!
    nullify ( from%value4, from%mask4 ) ! Don't deallocate during destroy!
    call destroyVectorQuantityValue ( from, destroyMask=.true. )
  end subroutine MoveVectorQuantity

  ! --------------------------------------------  MultiplyVectors  -----
  subroutine MultiplyVectors ( X, Y, Z, Quant, Inst )
  ! If Z is present, destroy Z and clone a new one from X, then
  ! Z = X # Y where # means "element-by-element"; otherwise X = X # Y

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    type(Vector_T), intent(in) :: Y
    type(Vector_T), intent(inout), optional, target :: Z
    integer, intent(in), optional :: Quant, Inst  ! If Quant is present,
    !  only that quantity is multiplied.  If furthermore Inst is present,
    !  only that instance is multiplied.  If Inst is present but Quant
    !  is not, the entire vector is multiplied.
    ! Local Variables:
    integer :: I                        ! Subscript and loop inductor
    type(Vector_T), pointer :: Result   ! associated to either X or Z
    ! Executable statements:
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
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

  ! ---------------------------------------- NullifyVectorTemplate -----
  subroutine NullifyVectorTemplate ( V )
    ! Given a vector template, nullify all the pointers associated with it
    type ( VectorTemplate_T ), intent(out) :: V

    ! Executable code
    v%name = 0
    nullify ( v%quantities )
  end subroutine NullifyVectorTemplate

  ! ------------------------------------------- NullifyVectorValue -----
  subroutine NullifyVectorValue ( V )
    ! Given a vector value, nullify all the pointers associated with it
    type ( VectorValue_T ), intent(out) :: V

    ! Executable code
    call nullifyQuantityTemplate ( v%template )
    nullify ( v%values, v%value1, v%value3 )
    nullify ( v%mask, v%mask1, v%mask3 )
  end subroutine NullifyVectorValue

  ! ---------------------------------------------- NullifyVector -----
  subroutine NullifyVector ( V )
    ! Given a vector, nullify all the pointers associated with it
    type ( Vector_T ), intent(out) :: V

    ! Executable code
    v%name = 0
    v%globalUnit = PHYQ_Invalid
    call nullifyVectorTemplate ( v%template )
    nullify ( v%quantities )
  end subroutine NullifyVector

  ! --------------------------------------------------  PowVector  -----
  subroutine PowVector ( X, POWER )
  ! Y = A / X if Y is present, else X = A / X.

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    real(rv), intent(in) :: POWER
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    do i = 1, size(x%quantities)
      x%quantities(i)%values = x%quantities(i)%values ** power
    end do
  end subroutine PowVector

  ! ------------------------------------------  ReciprocateVector  -----
  subroutine ReciprocateVector ( X, A, Y )
  ! Y = A / X if Y is present, else X = A / X.

    ! Dummy arguments:
    type(Vector_T), intent(inout), target :: X
    real(rv), intent(in) :: A
    type(Vector_T), intent(inout), optional, target :: Y
    ! Local Variables:
    integer :: I              ! Subscript and loop inductor
    type(Vector_T), pointer :: Z
    ! Executable statements:
    z => x
    if ( present(y) ) then
      if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Reciprocated vector has different template from original' )
      z => y
    end if
    do i = 1, size(x%quantities)
      z%quantities(i)%values = a / x%quantities(i)%values
    end do
  end subroutine ReciprocateVector

  ! --------------------------------------------  RemapVectorMask  -----
  subroutine RemapVectorMask ( Value )
    use Pointer_Rank_Remapping, only: Remap
    type ( vectorValue_t ) :: Value
    call remap ( value%mask1, value%mask, &
      & (/ value%template%noChans *     &
      &    value%template%noSurfs,      &
      &    value%template%noInstances * &
      &    value%template%noCrossTrack /) )
    call remap ( value%mask1, value%mask3, &
      & (/ value%template%noChans,      &
      &    value%template%noSurfs,      &
      &    value%template%noInstances * &
      &    value%template%noCrossTrack /) )
    call remap ( value%mask1, value%mask4, &
      & (/ value%template%noChans,      &
      &    value%template%noSurfs,      &
      &    value%template%noInstances,  &
      &    value%template%noCrossTrack /) )
  end subroutine RemapVectorMask

  ! -------------------------------------------  RemapVectorValue  -----
  subroutine RemapVectorValue ( Value )
    use Pointer_Rank_Remapping, only: Remap
    type ( vectorValue_t ) :: Value
    call remap ( value%value1, value%values, &
      & (/ value%template%noChans *     &
      &    value%template%noSurfs,      &
      &    value%template%noInstances * &
      &    value%template%noCrossTrack /) )
    call remap ( value%value1, value%value3, &
      & (/ value%template%noChans,      &
      &    value%template%noSurfs,      &
      &    value%template%noInstances * &
      &    value%template%noCrossTrack /) )
    call remap ( value%value1, value%value4, &
      & (/ value%template%noChans,      &
      &    value%template%noSurfs,      &
      &    value%template%noInstances,  &
      &    value%template%noCrossTrack /) )
  end subroutine RemapVectorValue

  ! -----------------------------------------  ReshapeVectorValue  -----
  ! Reshape source values to fit destination loosely
  ! destination must already be allocated
  ! Source values will come from either of
  ! (a) source%values, if present; or
  ! (b) an array of sourceValues
  ! If neither is present, we silently return an unchanged destination
  subroutine ReshapeVectorValue ( DESTINATION, SOURCE, SOURCEVALUES )
    type ( VectorValue_T), intent(out)                      :: DESTINATION
    type ( VectorValue_T), optional, pointer, intent(in)    :: SOURCE
    real(rv), dimension(:,:), optional, pointer, intent(in) :: SOURCEVALUES
    ! Internal variables
    integer :: i, j, k, n
    ! Executable
    if ( .not. associated(destination%values) ) return
    if ( present(source) ) then
      n = min( size(source%values), size(destination%values) )
      destination%value1(:n) = source%value1(:n)
    else if ( present(sourcevalues) ) then
      n = min( size(sourceValues), size(destination%values) )
      k = 0
      do j = 1, size(sourcevalues, 2)
        do i = 1, size(sourcevalues, 1)
          k = k + 1
          if ( k > n ) exit
          destination%value1(k) = sourcevalues(i, j)
        end do
      end do
    end if
  end subroutine ReshapeVectorValue

  ! ------------------------------------------------  ReverseMask  -----
  subroutine ReverseMask ( MASK, TO_Reverse, WHAT )
  ! Reverse bits of MASK indexed by elements of TO_Reverse.  Numbering of mask
  ! elements starts at one, not zero!  If TO_Reverse is absent, Reverse all of
  ! the bits of MASK.  If WHAT is absent, Reverse all bits.  If WHAT is
  ! present, Reverse only bits of MASK that correspond to "one" bits of WHAT.
  !
  ! In case it is not obvious, reverse turns "1" to "0", and "0" to "1"
  ! Effect: say the linAlg bits were set for elements {1,3,8,9} of a 10-element
  ! MASK. When we reverseMask the linAlg bits will be set for {2,4,5,6,7,10}
  ! All other bits and all other elements will be left unset
    character, intent(inout), dimension(:) :: MASK
    integer, intent(in), dimension(:), optional :: TO_Reverse
    integer, intent(in), optional :: WHAT
    integer :: MyWhat
    integer :: test
    logical, parameter :: DeeBug = .false.
    MyWhat = M_linAlg
    if ( present(what) ) myWhat = what
    if ( DeeBug ) then
      call output( 'Illustration: ' )
      test = 0
      call output (' effect on 0 is ')
      call output( ieor(test, myWhat ) )
      call blanks( 3 )
      test = 1
      call output (' effect on 1 is ')
      call output( ieor(test, myWhat ) )
    end if
    if ( present(to_Reverse) ) then
      mask(to_Reverse) = char(ieor(ichar(mask(to_Reverse)), myWhat))
    else
      mask = char(ieor(ichar(mask), myWhat))
    end if
  end subroutine ReverseMask

  ! ----------------------------------------  RmVectorFromDatabase  -----
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

  ! ------------------------------------------------  ScaleVector  -----
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
      if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, 'Scaled vector has different template from original' )
      z => y
    end if
    do i = 1, size(x%quantities)
      z%quantities(i)%values = a * x%quantities(i)%values
    end do
  end subroutine ScaleVector

  ! ----------------------------------------------------  SetMask  -----
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

  ! -----------------------------------------  SubtractFromVector  -----
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
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
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

  ! --------------------------------------------  SubtractVectors  -----
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
    if ( x%template%name /= y%template%name ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot subtract vectors having different templates" )
    call nullifyVector ( z ) ! for Sun's still useless compiler
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
    & noInstances, noSurfs, quantityType, molecule, signal, sideband, sayWhyNot )

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
    integer, optional, dimension(:), intent(IN) :: SIGNAL
    integer, optional, dimension(:), intent(IN) :: SIDEBAND
    logical, optional, intent(IN)               :: SAYWHYNOT

    ! Result
    logical :: ValidateVectorQuantity

    ! Executable code
    logical :: mySayWhyNot

    mySayWhyNot = .false.
    if ( present(sayWhyNot)) mySayWhyNot = sayWhyNot

    ValidateVectorQuantity = .true.

    if ( present(coherent) ) then
      if ( quantity%template%coherent .neqv. coherent ) then
        ValidateVectorQuantity=.FALSE.
        if ( mySayWhyNot ) then
          call output('Coherent quantity checked with incoherent', advance='yes')
          call output('quantity coherent? ', advance='no')
          call output(quantity%template%coherent, advance='yes')
          call output('check coherent? ', advance='no')
          call output(coherent, advance='yes')
        end if
        return
      end if
    end if

    if ( present(stacked) ) then
      if ( quantity%template%stacked .neqv. stacked ) then
        ValidateVectorQuantity=.FALSE.
        if ( mySayWhyNot ) then
          call output('stacked quantity checked with unstacked', advance='yes')
          call output('quantity stacked? ', advance='no')
          call output(quantity%template%stacked, advance='yes')
          call output('check stacked? ', advance='no')
          call output(stacked, advance='yes')
        end if
        return
      end if
    end if

    if ( present(regular) ) then
      if ( quantity%template%regular .neqv. regular ) then
        ValidateVectorQuantity=.FALSE.
        if ( mySayWhyNot ) then
          call output('Regular quantity checked with irregular', advance='yes')
          call output('quantity regular? ', advance='no')
          call output(quantity%template%regular, advance='yes')
          call output('check regular? ', advance='no')
          call output(regular, advance='yes')
        end if
        return
      end if
    end if

    if ( present(minorFrame) ) then
      if ( quantity%template%minorFrame .neqv. minorFrame ) then
        ValidateVectorQuantity=.FALSE.
        if ( mySayWhyNot ) then
          call output('Minor frame quantity checked with not', advance='yes')
          call output('quantity minor frame? ', advance='no')
          call output(quantity%template%minorFrame, advance='yes')
          call output('check minorFrame? ', advance='no')
          call output(minorFrame, advance='yes')
        end if
        return
      end if
    end if

    if ( present(majorFrame) ) then
      if ( quantity%template%majorFrame .neqv. majorFrame ) then
        ValidateVectorQuantity=.FALSE.
        if ( mySayWhyNot ) then
          call output('Major frame quantity checked with not', advance='yes')
          call output('quantity major frame? ', advance='no')
          call output(quantity%template%majorFrame, advance='yes')
          call output('check majorFrame? ', advance='no')
          call output(majorFrame, advance='yes')
        end if
        return
      end if
    end if

    if ( present(sideband) ) then
      ValidateVectorQuantity = any(quantity%template%sideband == sideband)
      if ( mySayWhyNot .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with different sideband', advance='yes')
        call output('quantity sideband ', advance='no')
        call output(quantity%template%sideband, advance='yes')
        call output('check sideband ', advance='no')
        call output(sideband, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(signal) ) then
      ValidateVectorQuantity = any(quantity%template%signal == signal)
      if ( mySayWhyNot .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with different signal', advance='yes')
        call output('quantity signal ', advance='no')
        call output(quantity%template%signal, advance='yes')
        call output('check signal ', advance='no')
        call output(signal, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(verticalCoordinate) ) then
      ValidateVectorQuantity=any(quantity%template%verticalCoordinate == verticalCoordinate)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with dif vert coord', advance='yes')
        call output('quantity vert coord ', advance='no')
        call output(quantity%template%verticalCoordinate, advance='yes')
        call output('check vert coord ', advance='no')
        call output(verticalCoordinate, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(frequencyCoordinate) ) then
      ValidateVectorQuantity=any(quantity%template%frequencyCoordinate == &
        & frequencyCoordinate)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with dif freq coord', advance='yes')
        call output('quantity freq coord ', advance='no')
        call output(quantity%template%frequencyCoordinate, advance='yes')
        call output('check freq coord ', advance='no')
        call output(frequencyCoordinate, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(noInstances) ) then
      ValidateVectorQuantity=any(quantity%template%noInstances == noInstances)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with dif num insts', advance='yes')
        call output('quantity num insts ', advance='no')
        call output(quantity%template%noInstances, advance='yes')
        call output('check noInstances ', advance='no')
        call output(noInstances, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(noSurfs) ) then
      ValidateVectorQuantity=any(quantity%template%noSurfs == noSurfs)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with dif num surfs', advance='yes')
        call output('quantity num surfs ', advance='no')
        call output(quantity%template%noInstances, advance='yes')
        call output('check noSurfs ', advance='no')
        call output(noSurfs, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(quantityType) ) then
      ValidateVectorQuantity=any(quantity%template%quantityType == quantityType)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with wrong type', advance='yes')
        call output('quantity type ', advance='no')
        call output(quantity%template%quantityType, advance='yes')
        call output('check quantityType ', advance='no')
        call output(quantityType, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

    if ( present(molecule) ) then
      ValidateVectorQuantity=any(quantity%template%molecule == molecule)
      if ( (mySayWhyNot) .and. .not. ValidateVectorQuantity ) then
        call output('quantity checked with wrong molecule', advance='yes')
        call output('quantity molecule ', advance='no')
        call output(quantity%template%molecule, advance='yes')
        call output('check molecule ', advance='no')
        call output(molecule, advance='yes')
      end if
      if ( .not. ValidateVectorQuantity) return
    end if

  end function ValidateVectorQuantity

  ! ------------------------------------------ VectorMemoryInUse -------
  integer function VectorMemoryInUse ( V )
  ! Report the total number of real(rv) elements in all quantities in vector V.
    type (Vector_T), intent(in) :: V
    integer :: I

    vectorMemoryInUse = 0
    do i = 1, size(v%quantities)
      if ( associated(v%quantities(i)%value1) ) &
        & vectorMemoryInUse = vectorMemoryInUse + size(v%quantities(i)%value1)
    end do

  end function VectorMemoryInUse

  ! ----------------------------------------- VectorsMemoryInUse -------
  integer function VectorsMemoryInUse ( D )
  ! Report the total number of real(rv) elements in all vectors in the
  ! database.
    type (Vector_T),  dimension(:) :: D
    integer :: I

    vectorsMemoryInUse = 0
    do i = 1, size(d)
      vectorsMemoryInUse = vectorsMemoryInUse + vectorMemoryInUse ( d(i) )
    end do

  end function VectorsMemoryInUse

! =====     Private Procedures     =====================================
  subroutine CreateValues ( Vector, highBound, lowBound )
  ! Allocate space for the values of a vector.
    type(Vector_T), intent(inout) :: Vector
    logical, intent(in), optional :: HIGHBOUND
    logical, intent(in), optional :: LOWBOUND
    integer :: QTY
    logical :: MYHIGHBOUND, MYLOWBOUND
    real(rv), parameter :: MYHUGE = 1.0e15
    character(63) :: What1, What2
    if ( vector%name == 0 ) then
      What1 = "Vector%"
    else
      call get_string ( vector%name, what1 )
      what1 = trim(what1) // "."
    end if
    myHighBound = .false.
    myLowBound = .false.
    if ( present ( highBound ) ) myHighBound = highBound
    if ( present ( lowBound ) ) myLowBound = lowBound
    do qty = 1, size(vector%quantities)
      if ( vector%quantities(qty)%template%name == 0 ) then
        what2 = "quantities(qty)"
      else
        call get_string ( vector%quantities(qty)%template%name, what2 )
      end if
      call createVectorValue ( vector%quantities(qty), trim(what1) // trim(what2) )
      if ( myHighBound ) then
        vector%quantities(qty)%values = myHuge
      else if ( myLowBound ) then
        vector%quantities(qty)%values = - myHuge
      else
        vector%quantities(qty)%values = 0.0_rv
      end if
    end do
  end subroutine

  subroutine get_string ( STRING, STRING_TEXT, CAP, STRIP, NOERROR, IERR, &
    & START, END )
  ! because get_string has the pernicious habit of bombing 
  ! if presented with 0 as its arg
  ! Args
    integer, intent(in) :: STRING
    character(len=*), intent(out) :: STRING_TEXT
    logical, intent(in), optional :: CAP
    logical, intent(in), optional :: STRIP
    logical, intent(in), optional :: NOERROR
    integer, intent(out), optional :: IERR
    integer, intent(in), optional :: START
    integer, intent(in), optional :: END
  string_text = 'undefined'
  if ( isStringInTable(string) ) &
    & call get_String_Rude( STRING, STRING_TEXT, CAP, STRIP, NOERROR, IERR, &
    & START, END )
  end subroutine get_string

!=======================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module VectorsModule
!=======================================================================

!
! $Log$
! Revision 2.211  2021/06/10 23:43:11  pwagner
! Added BinNumber and MAF components to vector qties
!
! Revision 2.210  2018/05/11 21:26:27  pwagner
! Moved M_ mask bit fields to MLSCommon
!
! Revision 2.209  2018/02/27 00:50:33  livesey
! Added the supportedInstrumentModule argument to the various search routines to support A-SMLS
!
! Revision 2.208  2017/11/03 19:57:59  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.207  2016/11/15 19:29:55  pwagner
! Prevent error from 0-size arrays when remapping pointer rank
!
! Revision 2.206  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.205  2016/05/27 00:14:05  vsnyder
! Publish RV because this seems like a logical place to get it
!
! Revision 2.204  2016/04/13 00:47:18  vsnyder
! Make AreEqual generic.  Add AreUnEqual, operator(==), and operator(/=)
!
! Revision 2.203  2015/09/24 18:46:09  pwagner
! Prevent an obscure cause of crashing when index too big in GetVectorQtyByTemplateIndex
!
! Revision 2.202  2015/09/17 22:55:23  pwagner
! Will now complain, quit in GetVectorQtyByTemplateIndex if index too big
!
! Revision 2.201  2015/07/14 23:20:47  pwagner
! Prevent more crashes due to get_string
!
! Revision 2.200  2015/06/19 00:35:49  pwagner
! Setup quick changes to print when destroying vectors
!
! Revision 2.199  2015/06/04 01:57:16  vsnyder
! Add contiguous attribute, avoid using remapped values
!
! Revision 2.198  2015/05/01 02:09:28  vsnyder
! Spiff a dump
!
! Revision 2.197  2015/04/29 00:53:56  vsnyder
! Spiff the dump
!
! Revision 2.196  2015/03/28 01:47:29  vsnyder
! Added 4-d Values and Mask.  Additional dimension extent is NoCrossTrack.
! Added stuff to trace allocate/deallocate addresses.
!
! Revision 2.195  2014/11/08 01:05:15  pwagner
! CreateVectorValue relies on DestroyVectorQuantityValue to deallocate its values
!
! Revision 2.194  2014/10/02 22:08:35  vsnyder
! Default initialize all components of Vector*_T
!
! Revision 2.193  2014/09/05 00:20:54  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Remove some
! debugging cruft that shouldn't have been checked in.
!
! Revision 2.192  2014/09/04 20:23:41  pwagner
! Turn off all those debugging print statements
!
! Revision 2.191  2014/08/19 00:29:59  vsnyder
! Add VectorMemoryInUse, VectorsMemoryInUse
!
! Revision 2.190  2014/07/18 23:13:24  pwagner
! Aimed for consistency in names passed to allocate_test; added allocationName field to datatype
!
! Revision 2.189  2014/05/20 23:53:19  vsnyder
! Add the POINTER attribute to the SourceValues argument to ReshapeVectorValue
! to work around bugs in compilers that have not properly implemented the
! feature that disassociated pointer actual arguments corresponding to
! optional non-pointer dummy arguments are absent.
!
! Revision 2.188  2014/02/01 00:16:43  pwagner
! Don't dump values for every quantity when failed to find one in GetVectorQuantityIndexByType
!
! Revision 2.187  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.186  2013/10/25 23:05:42  pwagner
! Fixed error that caused dumping hyperslab to segment fault if rank /= 2
!
! Revision 2.185  2013/09/27 00:41:14  pwagner
! Fixed bug in GatherVectorQuantity
!
! Revision 2.184  2013/09/25 00:58:15  pwagner
! Added a gather operation for Vector quantities
!
! Revision 2.183  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.182  2013/08/17 00:17:24  pwagner
! Guard against certain crashes with non-Aura datasets
!
! Revision 2.181  2013/08/16 02:28:38  vsnyder
! Don't reallocate VALUES if it's the right shape
!
! Revision 2.180  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.179  2013/07/13 00:00:15  vsnyder
! Define assignment, producing error message, for vector value
!
! Revision 2.178  2013/06/12 02:14:37  vsnyder
! Cruft removal
!
! Revision 2.177  2013/05/22 20:24:21  pwagner
! Dump procedures accept and obey options arg
!
! Revision 2.176  2013/05/16 18:16:58  pwagner
! Made GetVectorQuantityIndexByName generic--can supply character-valued arg
!
! Revision 2.175  2013/03/15 00:01:31  pwagner
! Added function AreEqual for vectors
!
! Revision 2.174  2013/02/21 21:21:42  pwagner
! Added options arg when dumping quantity or vector masks
!
! Revision 2.173  2013/02/01 23:40:55  vsnyder
! Add Where argument to CreateVectorValue.  Nullify correct pointers in
! DestroyVectorValue.
!
! Revision 2.172  2012/12/04 00:11:55  pwagner
! Improved comments
!
! Revision 2.171  2012/11/09 00:59:08  pwagner
! Added ReshapeVectorValue
!
! Revision 2.170  2012/10/30 22:07:28  pwagner
! Now able to clear full 2d quantity mask
!
! Revision 2.169  2012/10/29 17:42:31  pwagner
! Added optional args to CloneVectorQuantity, DestroyVectorQuantityValue
!
! Revision 2.168  2012/10/11 21:01:02  pwagner
! Print quantityName instead of moduleName during Dump
!
! Revision 2.167  2012/07/31 00:33:40  vsnyder
! Use Test_Allocate instead of explicit testing followed by MLSMessage.
! Add ForWhom argument to DestroyVectorQuantityMask and
! DestroyVectorQuantityValue.
!
! Revision 2.166  2012/07/10 03:55:49  vsnyder
! Improve comments about VALUE
!
! Revision 2.165  2012/07/07 02:40:47  vsnyder
! Add DestroyMask argument to DestroyVectorQuantityValue
!
! Revision 2.164  2012/07/07 02:02:34  vsnyder
! Add MASK1 and MASK3.  Make MASK and MASK3 rank remappings of MASK1.
! Add VALUE1 and VALUE3.  Make VALUES and VALUE3 rank remappings of VALUE1.
! Add low-level abstractions for creating and destroying masks and values
! for a single vector quantity.
!
! Revision 2.163  2012/05/24 20:32:56  vsnyder
! Change details level for dumping vector quantity templates
!
! Revision 2.162  2012/04/20 01:26:38  vsnyder
! Add DotVectorsMaybeMasked, norms, dump quantity norms
!
! Revision 2.161  2012/03/28 00:55:22  vsnyder
! Indicate mask is dumped in hex
!
! Revision 2.160  2012/02/23 00:08:55  vsnyder
! Don't dump molecule names if quantity type is not vmr
!
! Revision 2.159  2012/02/13 23:21:47  pwagner
! Print moleccule when dumping quantity
!
! Revision 2.158  2012/01/09 22:30:42  pwagner
! More info about dumped vector template
!
! Revision 2.157  2011/12/17 00:37:32  vsnyder
! Add MoveVectorQuantity
!
! Revision 2.156  2011/11/11 00:32:29  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.155  2011/11/01 22:55:48  honghanh
! Remove the nullify statement of vector%template%quantities
! from DestroyVectorInfo
!
! Revision 2.154  2011/10/25 18:08:07  pwagner
! Capitalize USEd items
!
! Revision 2.153  2011/08/02 16:51:40  honghanh
! Re-indent 2 lines of code to make indentation correct.
!
! Revision 2.152  2011/03/15 22:49:35  pwagner
! Added reverseMask; changed default dumpMask behavior to use DUMP_2D_INTEGER
!
! Revision 2.151  2011/03/02 02:15:25  vsnyder
! Make QuantityTemplate_t public, for F95 compatibility
!
! Revision 2.150  2010/05/24 14:48:04  honghanh
! Add comment to GetVectorQtyByTemplateIndex
!
! Revision 2.149  2010/04/28 00:12:20  pwagner
! Correct maskBitNames used in dumping mask
!
! Revision 2.148  2010/04/22 23:38:57  pwagner
! Added new Ignore masking bit
!
! Revision 2.147  2010/02/25 18:07:14  pwagner
! Added extra dump when about to bomb
!
! Revision 2.146  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.145  2009/10/27 23:28:46  pwagner
! Generic Diff must be public, too, at least NAG thinks so
!
! Revision 2.144  2009/10/27 22:16:38  pwagner
! New api for dump vector quantity--drops 'clean', adds 'options'
!
! Revision 2.143  2009/10/26 17:08:44  pwagner
! Added DiffVectorQuantities
!
! Revision 2.142  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.141  2009/06/16 17:17:55  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.140  2009/05/08 00:39:07  pwagner
! Shows MaskBitNames when dumping Mask bits
!
! Revision 2.139  2009/04/30 22:14:25  pwagner
! name of bit in MaskVectorQty and isVectorQtyMasked now mandatory
!
! Revision 2.138  2009/01/16 23:30:58  vsnyder
! Spiff up a dump
!
! Revision 2.137  2008/12/17 02:57:28  vsnyder
! Dump template with vector quantity if details gt 1
!
! Revision 2.136  2008/11/24 19:36:57  pwagner
! Improved comments, dumps; removed unused variables
!
! Revision 2.135  2008/11/06 21:51:08  pwagner
! Fill method swapValues swaps values between two quantities
!
! Revision 2.134  2008/08/27 19:58:30  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.133  2008/06/09 20:33:59  vsnyder
! Repair some broken comments
!
! Revision 2.132  2008/06/05 02:06:06  vsnyder
! Comments about Aux grids
!
! Revision 2.131  2007/10/09 16:45:02  pwagner
! Corrected declaration of myDetails
!
! Revision 2.130  2007/10/09 00:29:42  pwagner
! Added optional DETAILS arg to some dumps
!
! Revision 2.129  2007/10/06 00:01:28  vsnyder
! Remove unnecessary target attribute
!
! Revision 2.128  2007/10/03 20:51:18  vsnyder
! Add CheckVectorQuantityForNaN
!
! Revision 2.127  2007/04/03 17:41:58  vsnyder
! Revise how allocation status is tested.  Reallocate VectorsDatabase with
! zero size after destroying it.
!
! Revision 2.126  2006/08/05 02:11:58  vsnyder
! Add ForWhom argument to ConstructVectorTemplate
!
! Revision 2.125  2006/08/03 01:10:06  vsnyder
! Put l2cf names in leak track database
!
! Revision 2.124  2006/07/27 03:55:56  vsnyder
! Print summaries if negative details levels, for leak detection
!
! Revision 2.123  2006/06/06 18:54:48  vsnyder
! Spiff up a dump
!
! Revision 2.122  2006/05/23 21:43:34  vsnyder
! Add CLEAR option to some dumps
!
! Revision 2.121  2006/03/22 02:16:28  vsnyder
! Add Vector argument to DumpVectorQuantity just to get its name
!
! Revision 2.120  2006/02/23 00:55:28  vsnyder
! Add NoErr optional argument to GetVectorQuantityIndexByName
!
! Revision 2.119  2006/01/21 00:03:12  livesey
! Added DumpNiceMaskSummary
!
! Revision 2.118  2005/09/02 20:33:22  vsnyder
! Correct error messages in GetVectorQuantityIndexByType
!
! Revision 2.117  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.116  2005/03/03 02:14:01  vsnyder
! Spiff up some dumps
!
! Revision 2.115  2005/01/07 00:36:51  vsnyder
! Remove unused declarations
!
! Revision 2.114  2004/10/19 02:21:20  livesey
! Logical flaw in radiometer/signal vector querying
!
! Revision 2.113  2004/10/07 23:12:19  vsnyder
! Polish up Dump_Vector_Value for use in ForwardModelVectorTools
!
! Revision 2.112  2004/06/16 22:31:28  vsnyder
! Account for mask in DivideVectors, exchange order of first two arguments
!
! Revision 2.111  2004/06/16 01:18:39  vsnyder
! Add DivideVectors, incomplete comments in TOC about other routines
!
! Revision 2.110  2004/06/10 00:57:47  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.109  2004/05/01 04:07:44  vsnyder
! Rearranged some dumping stuff
!
! Revision 2.108  2004/01/30 23:27:59  livesey
! Added ReciprocateVector, PowVector and an optional argument to
! ClearVector
!
! Revision 2.107  2004/01/24 01:03:04  livesey
! Added allowNameMismatch argument to CopyVector
!
! Revision 2.106  2004/01/23 05:36:51  livesey
! Added DoVectors/QuantitiesMatch
!
! Revision 2.105  2003/09/15 23:28:50  vsnyder
! Remove unused private module variable
!
! Revision 2.104  2003/09/15 17:45:37  livesey
! Added target declaration for fussy intel compiler
!
! Revision 2.103  2003/08/27 20:06:42  livesey
! Bug fix in MaskVectorQty
!
! Revision 2.102  2003/06/20 19:33:53  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.101  2003/06/03 20:47:05  livesey
! Typo bug fix
!
! Revision 2.100  2003/06/03 19:23:03  livesey
! Added check to see that vector has not been destroyed
!
! Revision 2.99  2003/05/29 16:36:29  livesey
! New reflector argument to some of the GetVectorQuantity....
!
! Revision 2.98  2003/05/13 04:47:18  livesey
! Added noValues argument to CreateVector
!
! Revision 2.97  2003/05/12 02:05:27  livesey
! Added InflateVectorTemplateDatabase and InflateVectorDatabase
!
! Revision 2.96  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.95  2003/04/04 22:54:58  livesey
! New mask bits
!
! Revision 2.94.2.3  2003/04/15 23:13:38  vsnyder
! Pass 'clean' option through dump_vector_quantity
!
! Revision 2.94.2.2  2003/03/07 23:51:17  vsnyder
! Copy the quantity index in CopyVector
!
! Revision 2.94.2.1  2003/03/06 23:25:16  vsnyder
! Add Dump_Vector_Quantity
!
! Revision 2.94  2002/11/22 12:57:09  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.93  2002/10/19 18:53:26  livesey
! Changed from huge to our own value (temporarily?)
!
! Revision 2.92  2002/10/17 18:18:25  livesey
! Added low/high bound stuff
!
! Revision 2.91  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.90  2002/09/26 18:01:08  livesey
! Made GetVectorQuantity... more forgiving in the case of l_vmr (can have
! radiometer wrong if molecule is extinction).
!
! Revision 2.89  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.88  2002/09/11 14:06:12  livesey
! Bug fix in CopyVector
!
! Revision 2.87  2002/08/08 22:06:33  vsnyder
! Add M_Tikhonov
!
! Revision 2.86  2002/08/04 15:55:56  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.85  2002/07/22 03:26:15  livesey
! Added CheckIntegrity
!
! Revision 2.84  2002/07/01 23:51:30  vsnyder
! Plug a memory leak
!
! Revision 2.83  2002/05/17 17:56:01  livesey
! More checks in ValidateVectorQuantity
!
! Revision 2.82  2002/05/14 00:28:04  livesey
! More informative messages in GetVectorQuantityByType
!
! Revision 2.81  2002/04/22 20:54:29  vsnyder
! Add a 'scale' argument to AddToVector
!
! Revision 2.80  2002/03/13 22:00:16  livesey
! Changed m_explicitFill to m_fill
!
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

