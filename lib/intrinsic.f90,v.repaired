! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module intrinsic

! Intrinsic constants needed by Init_Tables_Module, DeclarationTable, etc.

! Declaring the definitions is handled by the tree walker.

  implicit none
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! A "spec_def" vertex may be decorated with (sums of) the following flags:
  integer, parameter :: NO_DUP = 1         ! Duplicate fields prohibited
  integer, parameter :: ALL_FIELDS = 2     ! All fields required
  integer, parameter :: NO_POSITIONAL = 4  ! Positional fields prohibited
! A "field_type", "field_spec" or "dot" vertex may be decorated with (sums of)
! the following flags:
  integer, parameter :: NO_ARRAY = 1       ! Field must be scalar
  integer, parameter :: ARRAY_ARRAY = 2    ! Array elements as arrays OK
  integer, parameter :: REQ_FLD = 4        ! Required field
  integer, parameter :: EMPTY_OK = 8       ! Field can have empty value
  integer, parameter :: EXPR_OK = 16       ! Expr of dots is OK
  integer, parameter :: U = 64             ! U*PHYQ_... requires those units
! A "section" vertex may be decorated with the following flag:
  integer, parameter :: NO_CHECK_EQ = 1    ! Don't check whether the section's
                                           ! A=B contents are allowed.

  integer, parameter :: T_First             = 1
! Types of entities in the declaration table:
  integer, parameter :: T_Unknown           = t_first
  integer, parameter :: T_Empty             = t_unknown + 1
  integer, parameter :: T_Do_label          = t_empty + 1 ! DO, CASE or IF
  integer, parameter :: T_Enum_name         = t_do_label + 1
  integer, parameter :: T_Field_name        = t_enum_name + 1
  integer, parameter :: T_Function_name     = t_field_name + 1
  integer, parameter :: T_Label             = t_function_name + 1 ! of a spec
  integer, parameter :: T_Param_name        = t_label + 1 
  integer, parameter :: T_Phys_unit_name    = t_param_name + 1    ! PHYS_....
  integer, parameter :: T_Section_name      = t_phys_unit_name + 1
  integer, parameter :: T_Spec_name         = t_section_name + 1
  integer, parameter :: T_Tree_name         = t_spec_name + 1     ! e.g. n_plus
  integer, parameter :: T_Type_name         = t_tree_name + 1
  integer, parameter :: T_Unit_name         = t_type_name + 1
  integer, parameter :: T_Variable_name     = t_unit_name + 1
! Result type from Vector_Qty_Expr
  integer, parameter :: T_A_DOT_B           = t_variable_name + 1
! Declaration table types of names of LHS entities in Algebra_m
  integer, parameter :: T_Exprn             = t_a_dot_b + 1 ! Scalar
  integer, parameter :: T_Exprn_m           = t_exprn + 1         ! Matrix
  integer, parameter :: T_Exprn_v           = t_exprn_m + 1       ! Vector
! Data types that don't have enumerated literals:
  integer, parameter :: T_Numeric           = t_exprn_v + 1
  integer, parameter :: T_Numeric_range     = t_numeric + 1
  integer, parameter :: T_String            = t_numeric_range + 1
  integer, parameter :: T_String_range      = t_string + 1
! Enumeration types:
  integer, parameter :: T_Boolean           = t_string_range + 1
  integer, parameter :: T_Instrument        = t_boolean + 1
  integer, parameter :: T_Polarization      = t_instrument + 1
  integer, parameter :: Last_Intrinsic_Type = t_polarization

! We don't define any fields here, but here's the first index:
  integer, parameter :: Field_First = 1

! Abstract physical quantities:
  integer, parameter :: First_PhyQ = 0
  integer, parameter :: PhyQ_Invalid =         first_phyq ! Invalid unit given by user
  integer, parameter :: PhyQ_Dimensionless =   phyq_invalid+1     ! Dimensionless quantity
  integer, parameter :: PhyQ_Length =          phyq_dimensionless+1  ! Default meters
  integer, parameter :: PhyQ_Time =            phyq_length+1         ! Default seconds
  integer, parameter :: PhyQ_Pressure =        phyq_time+1        !  Default millibars
  integer, parameter :: PhyQ_Temperature =     phyq_pressure+1    ! Default Kelvins
  integer, parameter :: PhyQ_Vmr =             phyq_temperature+1 ! Default parts-per-one
  integer, parameter :: PhyQ_Angle =           phyq_vmr+1         ! Default degrees
  integer, parameter :: PhyQ_Mafs =            phyq_angle+1       ! Default MAFs
  integer, parameter :: PhyQ_Mifs =            phyq_mafs+1        ! Default MIFs
  integer, parameter :: PhyQ_Frequency =       phyq_mifs+1        ! Default MHz
  integer, parameter :: PhyQ_Zeta =            phyq_frequency+1   ! log10(pressure/hPa)
  integer, parameter :: PhyQ_Velocity =        phyq_zeta+1        ! Default meters/second
  integer, parameter :: PhyQ_Extinction =      phyq_velocity+1    ! Default 1/meters
  integer, parameter :: PhyQ_Icedensity =      phyq_extinction+1  ! Default g/meters^3
  integer, parameter :: PhyQ_Colmabundance =   phyq_icedensity+1  ! Default log10 g/meters^3
  integer, parameter :: PhyQ_Pctrhi =          phyq_colmabundance+1 ! default %RHI
  integer, parameter :: PhyQ_Gauss =           phyq_pctrhi + 1
  integer, parameter :: PhyQ_Profiles =        phyq_gauss + 1
  integer, parameter :: Last_phyq = phyq_profiles
  integer :: PhyQ_Indices(first_phyq:last_phyq)

! Enumeration literals:
  integer, parameter :: FIRST_LIT       = 1
  integer, save      :: LAST_AUTO_LIT   = 0 ! INIT_TABLES_MODULE should reset 
! Don't edit the following file directly -- it is generated automatically
! based on the file lit_names.txt (which is the file you ought to edit).
  include 'lit_parm.f9h'

  ! Specifications
  integer, parameter :: Spec_First = 1
  integer, parameter :: S_TIME          = Spec_First
  integer, parameter :: Last_Intrinsic_Spec = S_Time

  ! The following parameters are for building trees:
  integer, parameter :: BEGIN = -1  ! Start of a tree
  integer, parameter :: D = 1000000 ! Decoration, D*1 on n_func_def requires
                                    ! uniform but unspecified argument types
  integer, parameter :: DU = d*u    ! DU*PHYQ_... on n_field_type requires those
                                    ! units.  DU*type on n_func_def specifies
                                    ! function result type, else the result
                                    ! type is the same as the first argument.
  integer, parameter :: F = 1000    ! Field index
  integer, parameter :: G = 2000    ! Function index
  integer, parameter :: L = 3000    ! Lit index
  integer, parameter :: N = 0       ! Tree index
  integer, parameter :: NADP = n+d*(all_fields+no_dup+no_positional)
  integer, parameter :: NC = n+d*no_check_eq
  integer, parameter :: ND = n+d*no_dup
  integer, parameter :: NDP = n+d*(no_dup+no_positional)
  integer, parameter :: NDR = n+d*(no_dup+req_fld)
  integer, parameter :: NP = n+d*no_positional
  integer, parameter :: NR = n+d*req_fld
  integer, parameter :: NRS = n+d*(no_array+req_fld)
  integer, parameter :: NS = n+d*no_array
  integer, parameter :: P = 4000    ! Parameter index
  integer, parameter :: S = 5000    ! Spec index
  integer, parameter :: T = 6000    ! Type index
  integer, parameter :: Z = 7000    ! Section index

  ! Tables used for type checking:
  integer, save, pointer, dimension(:) :: Data_Type_Indices  => NULL()
  integer, save, pointer, dimension(:) :: Field_Indices      => NULL()
  integer, save, pointer, dimension(:) :: Func_Indices       => NULL()
  integer, save, pointer, dimension(:) :: Lit_Indices        => NULL()
  integer, save, pointer, dimension(:) :: Parm_Indices       => NULL()
  integer, save, pointer, dimension(:) :: Section_Indices    => NULL()
  integer, save, pointer, dimension(:) :: Spec_Indices       => NULL()
  
  ! Private procedures
  private :: allocate_test, deallocate_test

contains ! =====     Public procedures     =============================
! -----------------------------------------------  INIT_INTRINSIC  -----
  subroutine INIT_INTRINSIC ( N_DATA_TYPE_INDICES, N_FIELD_INDICES, &
    & N_LIT_INDICES, FIRST_PARM_INDEX, LAST_PARM_INDEX, N_SECTION_INDICES, &
    & N_SPEC_INDICES, N_FUNC_INDICES )

    ! This really belongs in make_tree, but "make depends" can't see it there
    ! (because of the "include"):

    ! use Allocate_Deallocate, only: Allocate_Test
    use Tree, only: ! Build_tree, Push_pseudo_terminal
    use Tree_Types, only: N_Dt_Def

    integer, intent(in) :: N_DATA_TYPE_INDICES
    integer, intent(in) :: N_FIELD_INDICES
    integer, intent(in) :: N_LIT_INDICES
    integer, intent(in) :: FIRST_PARM_INDEX, LAST_PARM_INDEX
    integer, intent(in) :: N_SECTION_INDICES
    integer, intent(in) :: N_SPEC_INDICES
    integer, intent(in) :: N_FUNC_INDICES

    ! Allocate the string index tables for the various categories of
    ! names
    call allocate_test ( data_type_indices, n_data_type_indices, &
      & 'DATA_TYPE_INDICES', moduleName )
    call allocate_test ( field_indices, n_field_indices, &
      & 'FIELD_INDICES', moduleName )
    call allocate_test ( lit_indices, n_lit_indices, &
      & 'LIT_INDICES', moduleName )
    call allocate_test ( parm_indices, last_parm_index, &
      & 'PARM_INDICES', moduleName, lowBound=first_parm_index )
    call allocate_test ( section_indices, n_section_indices, &
      & 'SECTION_INDICES', moduleName )
    call allocate_test ( spec_indices, n_spec_indices, &
      & 'SPEC_INDICES', moduleName )
    call allocate_test ( func_indices, n_func_indices, &
      & 'FUNC_INDICES', moduleName )

    ! Put intrinsic predefined identifiers into the symbol table.

    ! Put declaration-table names into symbol table
    data_type_indices(t_unknown) =         add_ident ( 'unknown' )
    data_type_indices(t_empty) =           add_ident ( 'empty' )
    data_type_indices(t_do_label) =        add_ident ( 'do label' )
    data_type_indices(t_enum_name) =       add_ident ( 'enum name' )
    data_type_indices(t_field_name) =      add_ident ( 'field name' )
    data_type_indices(t_function_name) =   add_ident ( 'function name' )
    data_type_indices(t_label) =           add_ident ( 'label' )
    data_type_indices(t_param_name) =      add_ident ( 'param name' )
    data_type_indices(t_phys_unit_name) =  add_ident ( 'phys unit name' )
    data_type_indices(t_section_name) =    add_ident ( 'section name' )
    data_type_indices(t_spec_name) =       add_ident ( 'spec name' )
    data_type_indices(t_tree_name) =       add_ident ( 'tree name' )
    data_type_indices(t_type_name) =       add_ident ( 'type name' )
    data_type_indices(t_unit_name) =       add_ident ( 'unit name' )
    data_type_indices(t_variable_name) =   add_ident ( 'variable name' )
    ! Put result type from Vector_Qty_Expr into the symbol table
    data_type_indices(t_a_dot_b) =         add_ident ( 'a.b' )
    ! Put type names of declarations of Algebra_m LHS's into the symbol table
    data_type_indices(t_exprn) =           add_ident ( 'exprn' )
    data_type_indices(t_exprn_m) =         add_ident ( 'exprn_m' )
    data_type_indices(t_exprn_v) =         add_ident ( 'exprn_v' )
    ! Put intrinsic non-enumeration type names into the symbol table
    data_type_indices(t_numeric) =         add_ident ( 'numeric' )
    data_type_indices(t_numeric_range) =   add_ident ( 'numeric_range' )
    data_type_indices(t_string) =          add_ident ( 'string' )
    data_type_indices(t_string_range) =    add_ident ( 'string_range' )
    ! Put intrinsic enumeration type names into the symbol table
    data_type_indices(t_boolean) =         add_ident ( 'boolean' )
    data_type_indices(t_instrument) =      add_ident ( 'instrument' )
    data_type_indices(t_polarization) =    add_ident ( 'polarization' )
    ! Put intrinsic enumeration literals into the symbol table:
! Don't edit the following file directly--it is generated automatically
! based on the file lit_names.txt
    include 'lit_add.f9h'

    ! Put spec names into the symbol table
    spec_indices(s_time) =                 add_ident ( 'time' )

    ! Put abstract physical quantities into the symbol table
    phyq_indices(phyq_invalid) =           add_ident ( 'invalid' )
    phyq_indices(phyq_dimensionless) =     add_ident ( 'dimensionless' )
    phyq_indices(phyq_length) =            add_ident ( 'length' )
    phyq_indices(phyq_time) =              add_ident ( 'time' )
    phyq_indices(phyq_pressure) =          add_ident ( 'pressure' )
    phyq_indices(phyq_temperature) =       add_ident ( 'temperature' )
    phyq_indices(phyq_vmr) =               add_ident ( 'vmr' )
    phyq_indices(phyq_angle) =             add_ident ( 'angle' )
    phyq_indices(phyq_mafs) =              add_ident ( 'mafs' )
    phyq_indices(phyq_mifs) =              add_ident ( 'mifs' )
    phyq_indices(phyq_frequency) =         add_ident ( 'frequency' )
    phyq_indices(phyq_velocity) =          add_ident ( 'velocity' )
    phyq_indices(phyq_zeta) =              add_ident ( 'zeta' )
    phyq_indices(phyq_extinction) =        add_ident ( 'extinction' )
    phyq_indices(phyq_icedensity) =        add_ident ( 'iceDensity' )
    phyq_indices(phyq_colmabundance) =     add_ident ( 'colmabundance' )
    phyq_indices(phyq_pctrhi) =            add_ident ( 'pctrhi' )
    phyq_indices(phyq_gauss) =             add_ident ( 'gauss' )
    phyq_indices(phyq_profiles) =          add_ident ( 'profiles' )

  ! Definitions are represented by trees.  The notation in the comments
  ! for the trees is < root first_son ... last_son >.  This is sometimes
  ! called "Cambridge Polish Notation."  It was developed to represent
  ! LISP by McCarthy et. al. at MIT (in Cambridge, MA).

  ! Notice that in the argument for make_tree, the tree node id is at
  ! the END of the subtree, while in Cambridge Polish Notation it is at
  ! the BEGINNING of the subtree!

  ! Put the definition trees into the tree space before the parser runs.
  ! After the parsing is done, they're automatically "glued in" to the
  ! "left" of the trees that represent the input.  The tree-walker
  ! stumbles upon them in its normal course of operation, never really
  ! realizing they're special (because by then they're not).

  ! Start with the definitions of types. These are represented by trees of
  ! the form  < n_dt_def t_type_name l_lit ... l_lit >

    ! Define the intrinsic data types
    call make_tree ( (/ &
      begin, t+t_numeric, n+n_dt_def, &
      begin, t+t_numeric_range, n+n_dt_def, &
      begin, t+t_string, n+n_dt_def, &
      begin, t+t_string_range, n+n_dt_def /) )
    ! Define the enumerated types
    call make_tree ( (/ &
      begin, t+t_boolean, l+l_true, l+l_false, n+n_dt_def,   &
      begin, t+t_instrument, l+l_emls, l+l_umls, l+l_xptl1, n+n_dt_def, &
      begin, t+t_polarization, l+l_a, l+l_b, n+n_dt_def /) )

  contains
    ! ................................................  MAKE_TREE  .....
    include "make_tree.f9h"
    
  end subroutine INIT_INTRINSIC

  ! --------------------------------------------------  Add_Ident  -----
  integer function ADD_IDENT ( TEXT )
    use Symbol_Table, only: Enter_Terminal
    use Symbol_Types, only: T_Identifier
    character(len=*), intent(in) :: TEXT
    add_ident = enter_terminal ( text, t_identifier )
  end function ADD_IDENT

  ! -----------------------------------  DestroyTypeCheckerTables  -----
  subroutine DestroyTypeCheckerTables
    ! use Allocate_Deallocate, only: Deallocate_Test
    call deallocate_test ( data_type_indices, 'DATA_TYPE_INDICES', moduleName )
    call deallocate_test ( field_indices,     'FIELD_INDICES',     moduleName )
    call deallocate_test ( lit_indices,       'LIT_INDICES',       moduleName )
    call deallocate_test ( parm_indices,      'PARM_INDICES',      moduleName )
    call deallocate_test ( section_indices,   'SECTION_INDICES',   moduleName )
    call deallocate_test ( spec_indices,      'SPEC_INDICES',      moduleName )
  end subroutine

  ! ---------------------------------------------------  Get_Type  -----
  integer function Get_Type ( Type_Index )
    ! Return the string index for a type name.  The reason for the
    ! existence of this function is to pass it to procedures that otherwise
    ! do not need to access this module by use association.
    integer, intent(in) :: Type_Index
    if ( type_index >= lbound(data_type_indices,1) .and. &
       & type_index <= ubound(data_type_indices,1) ) then
      get_type = data_type_indices ( type_index )
    else
      get_type = 0
    end if
  end function Get_Type
  
  ! Try to avoid USE-ing Allocate_Deallocate to lift circular dependency
  subroutine deallocate_test ( indices, type_str, whereami )
    ! Args
    integer, pointer, dimension(:)          :: indices
    character(len=*), intent(in)            :: type_str, whereami
    ! Internal variables
    integer                                 :: status ! 0 success, ! 0 failure
    character(len=127)                      :: ermsg  ! any clue why?
    ! Executable
    deallocate( indices, stat=status, errmsg=ermsg )
    if ( status == 0 ) return
    print *, 'Failed to deallocate ' // type_str // ' in ' // whereami
    print *, trim(ermsg)
  end subroutine deallocate_test

  subroutine allocate_test ( indices, n, type_str, whereami, lowbound )
    ! Args
    integer, pointer, dimension(:)          :: indices
    integer, intent(in)                     :: n
    character(len=*), intent(in)            :: type_str, whereami
    integer, intent(in), optional           :: lowbound
    ! Internal variables
    integer                                 :: n1 ! low bound of indices
    integer                                 :: status ! 0 success, ! 0 failure
    character(len=127)                      :: ermsg  ! any clue why?
    ! Executable
    n1 = 1
    if ( present(lowbound) ) n1 = lowbound
    allocate( indices(n1:n), stat=status, errmsg=ermsg )
    if ( status == 0 ) return
    print *, 'Failed to allocate ' // type_str // ' in ' // whereami
    print *, trim(ermsg)
  end subroutine allocate_test

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module INTRINSIC

! $Log$
! Revision 2.76  2020/07/28 20:31:43  vsnyder
! Added ARRAY_ARRAY type-checker flag. Changed value of U flag
!
! Revision 2.75  2019/08/19 22:00:23  pwagner
! Avoid USE-ing Allocate_Deallocate due to circular dependency
!
! Revision 2.74  2016/10/21 23:28:20  vsnyder
! Remove unused USE name
!
! Revision 2.73  2014/05/20 22:16:24  vsnyder
! Don't go out of bounds in Get_Type
!
! Revision 2.72  2014/03/20 01:38:29  vsnyder
! Unify types in Intrinsic instead of having a separate system in
! Declaration_Table.
!
! Revision 2.71  2013/12/12 01:59:44  vsnyder
! Change type of debug from logical to integer, add 'unknown' type
!
! Revision 2.70  2013/10/09 01:05:18  vsnyder
! Comments to clarify use of D and DU
!
! Revision 2.69  2013/09/19 23:26:44  vsnyder
! Add Expr_OK
!
! Revision 2.68  2013/09/17 00:55:49  vsnyder
! Add A_Dot_B type
!
! Revision 2.67  2012/01/05 01:14:11  pwagner
! Added get_phyq function to phyq indices; also LAST_AUTO_LIT
!
! Revision 2.66  2011/04/18 19:28:05  vsnyder
! Add NC for unchecked parameters in sections
!
! Revision 2.65  2011/01/29 00:46:42  vsnyder
! Add units checking
!
! Revision 2.64  2010/02/04 23:03:25  vsnyder
! Remove PHYQ_LOGICEDENSITY
!
! Revision 2.63  2009/09/19 00:35:07  vsnyder
! Add phyq_logIceDensity
!
! Revision 2.62  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.61  2008/08/27 19:58:30  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.60  2006/03/23 01:50:12  vsnyder
! Add Empty_OK parameter
!
! Revision 2.59  2006/01/11 16:59:36  pwagner
! Abstract phys quant now colmabundance
!
! Revision 2.58  2005/12/29 01:09:49  vsnyder
! Add string for PHYQ_Profiles
!
! Revision 2.57  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.56  2005/04/19 19:13:10  livesey
! Changed mls1 to xptl1
!
! Revision 2.55  2004/11/17 20:23:09  vsnyder
! Add NRS and NS (scalar required) tags for fields
!
! Revision 2.54  2004/05/29 02:42:59  vsnyder
! Rearrange function definition stuff
!
! Revision 2.53  2004/01/14 18:50:25  vsnyder
! Stuff to support the Algebra section
!
! Revision 2.52  2004/01/09 07:24:59  livesey
! Added the fictitious instrument mls1
!
! Revision 2.51  2003/09/15 17:01:09  livesey
! Put the use statement back, it's needed because the includes need
! it and make doesn't spot it in there.
!
! Revision 2.49  2003/08/16 01:14:03  vsnyder
! Add optional 'polarization' field to 'radiometer' spec
!
! Revision 2.48  2003/07/08 00:16:08  livesey
! Added ndr
!
! Revision 2.47  2003/05/29 16:36:09  livesey
! Added ident for phyq_gauss
!
! Revision 2.46  2003/01/26 04:41:53  livesey
! Added phyq_profiles
!
! Revision 2.45  2003/01/07 23:43:59  livesey
! Added Gauss
!
! Revision 2.44  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.43  2002/04/10 17:42:59  pwagner
! Added pctrhi unit
!
! Revision 2.42  2001/10/04 22:12:57  pwagner
! Now includes files lit_add.f9h and lit_parm.f9h
!
! Revision 2.41  2001/10/03 18:32:44  vsnyder
! OOPS, forgot some lit_indices(...) = add_ident(...)
!
! Revision 2.40  2001/10/03 17:38:11  vsnyder
! OOPS, defined L_DNWT_SQ twice
!
! Revision 2.39  2001/10/03 17:36:47  vsnyder
! Add lits for DNWT quantities
!
! Revision 2.38  2001/10/02 23:39:39  vsnyder
! Add L_DegreesOfFreedom
!
! Revision 2.37  2001/09/17 23:14:14  livesey
! Bug fix, added name for t_instrument
!
! Revision 2.36  2001/09/17 22:53:23  livesey
! Added t_instrument, l_emls and l_umls
!
! Revision 2.35  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.34  2001/07/18 23:15:31  dwu
! rename l_radiusofearth as l_earthradius
!
! Revision 2.33  2001/07/17 18:53:31  jonathan
! remove earthradius (redundant quantity),jonathan/wu
!
! Revision 2.32  2001/07/13 19:04:56  dwu
! fix problem after adding lostransfunc
!
! Revision 2.31  2001/07/13 18:20:18  dwu
! add quantity losTransFunc
!
! Revision 2.30  2001/07/10 23:46:35  jonathan
! added l_icedensity, paul/jonathan
!
! Revision 2.29  2001/07/06 18:55:40  jonathan
! Modified for cloud model, Paul/Jonathan
!
! Revision 2.28  2001/05/31 22:07:33  livesey
! Updated cloud quantity types
!
! Revision 2.27  2001/05/31 20:27:37  livesey
! New vector type associated with cloud quantities.
!
! Revision 2.26  2001/05/29 22:45:45  livesey
! Moved some state vector component literals in from l2/init_tables_module.f90
!
! Revision 2.25  2001/05/10 23:26:45  livesey
! Added isotope ratio vector quantity type.
!
! Revision 2.24  2001/05/03 22:27:04  livesey
! Added l_heightoffset
!
! Revision 2.23  2001/04/26 16:18:39  livesey
! Fixed bug. All those nice integer arrays weren't initially nullified, whoops!
!
! Revision 2.22  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.21  2001/04/23 20:57:42  vsnyder
! Move the first spec (time) to 'intrinsic'
!
! Revision 2.20  2001/04/09 20:59:35  vsnyder
! Add C (for Celsius) unit and l_c name for it
!
! Revision 2.19  2001/04/09 14:56:04  perun
! Corrected enumeration literal definition for l_temperature
!
! Revision 2.18  2001/04/04 17:56:42  vsnyder
! Insert "USE TREE" because "make depends" can't see the one in "make_tree"
! (because of the "include").
!
! Revision 2.17  2001/04/04 17:21:12  pwagner
! Added extra use tree line to tweak dependencies
!
! Revision 2.16  2001/04/03 19:09:12  vsnyder
! Change the order of initialization to intrinsic, Molecules, MLSSignals.
! Use the revised make_tree.f9h, which requires revision of init...
! calling sequences.
!
! Revision 2.15  2001/03/17 02:23:40  livesey
! Bug fix, defined phyq_indices(phyq_velocity)
!
! Revision 2.14  2001/03/15 18:41:04  livesey
! Added some more, losvel etc.
!
! Revision 2.13  2001/03/15 07:37:35  livesey
! Added l_frequency
!
! Revision 2.12  2001/03/14 02:05:52  vsnyder
! Moved MLSSignals_m to mlspgs/lib.
!
! Revision 2.11  2001/03/06 22:41:59  livesey
! Minor changes
!
! Revision 2.10  2001/03/02 01:32:21  livesey
! Added some new PHYQs
!
! Revision 2.9  2001/02/22 23:57:38  vsnyder
! Remove ", public" from parameters, because default accessibility is public
!
! Revision 2.8  2001/02/09 18:37:37  vsnyder
! Add REQ_FLD flag for specification definitions
!
! Revision 2.7  2001/02/09 01:05:18  livesey
! Thought I'd done this
!
! Revision 2.6  2001/02/08 21:11:13  vsnyder
! Move "theta" from init_tables_module to intrinsic.
!
! Revision 2.5  2001/02/05 21:18:57  vsnyder
! Add parameters for type checking rules.
!
! Revision 2.4  2001/02/01 20:18:50  vsnyder
! Correct index and spelling for gph_precision and temperature_precision
!
! Revision 2.3  2001/02/01 01:23:18  vsnyder
! Account for the Molecules module
!
! Revision 2.2  2001/01/31 23:32:31  vsnyder
! Moved l_temperature l_temperature_prec l_ptan l_tangentheight l_sidebandratio
! l_scvel l_orbitinclination l_geodaltitude l_radiance l_scanresidual l_gph
! l_gph_precision l_refgph l_baseline l_extinction l_linewidth from L2's
! init_tables_module
!
! Revision 2.1  2000/10/11 18:24:39  vsnyder
! Initial entry
!
