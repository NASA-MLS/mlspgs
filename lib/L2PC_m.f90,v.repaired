! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module L2PC_m
  !=============================================================================

  ! This module contains data types etc. for dealing with the new EMLS L2PC
  ! files.  The first version dealt with ascii files, but later versions
  ! must be HDF5.

  use Allocate_deallocate, only: allocate_test, deallocate_test
  use Dump_0, only: dump
  use HessianModule_0, only: hessianElement_t, &
    & H_absent, h_sparse, h_full, h_unknown, &
    & Createblock, destroyblock
  use HessianModule_1, only: hessian_t, &
    & Copyhessianvalue, createblock, createemptyhessian, destroyhessian
  use Highoutput, only: outputnamedvalue
  use Intrinsic, only: l_adopted, l_channel, l_geodaltitude, l_none, l_vmr, &
    & L_radiance, l_none, l_intermediatefrequency, l_latitude, l_fieldazimuth, &
    & L_rows, l_columns, l_columnabundance, l_temperature, l_tscat, &
    & L_isotoperatio, l_calsidebandfraction, l_limbsidebandfraction, &
    & L_opticaldepth, l_elevoffset, &
    & Lit_indices, &
    & Phyq_colmabundance, phyq_dimensionless, phyq_pctrhi, &
    & Phyq_temperature, phyq_vmr
  use Manipulatevectorquantities, only: dovectorsmatch
  use MatrixModule_0, only: m_absent, m_banded, m_column_sparse, m_full, &
    & MatrixElement_t, m_unknown, destroyblock
  use MatrixModule_1, only: matrix_t, matrix_database_t, &
    & CopyMatrixValue, createBlock, createEmptyMatrix, &
    & DestroyMatrix, dump, dump_struct, &
    & FindBlock, getActualMatrixFromDatabase
  use MLSCommon, only: MLSFile_t
  use MLSFiles, only: dumpMLSFile => dump
  use MLSKinds, only: r8, r4
  use MLSMessagemodule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSFinds, only: findfirst
  use MLSSignals_m, only: getsignalname
  use MLSStringlists, only: optiondetail, switchdetail
  use MLSStrings, only: writeintstochars
  use Molecules, only: isextinction, l_rhi
  use Moretree, only: getlitindexfromstring, getstringindexfromstring
  use Output_m, only: newline, output
  use Parse_signal_m, only: parse_signal
  use Quantitytemplates, only: quantitytemplate_t, &
    & Addquantitytemplatetodatabase, copyquantitytemplate, &
    & Destroyquantitytemplatecontents, inflatequantitytemplatedatabase, &
    & Nullifyquantitytemplate, setupnewquantitytemplate
  use String_table, only: display_string, getString => get_string, &
    & isStringInTable
  use Toggles, only: switches
  use Tree, only: decoration, nsons, subtree
  use Vectorsmodule, only: vectortemplate_t, vector_t, &
    & Assignment(=), addvectortemplatetodatabase, &
    & Addvectortodatabase, constructvectortemplate, copyvector, createvector, &
    & Destroyvectorinfo, dump, nullifyvectortemplate

  implicit none
  private

  public :: AddBinSelectorToDatabase, AddL2PCToDatabase, AdoptVectorTemplate, &
    & binSelector_T, BinSelectors, CreateDefaultBinSelectors, &
    & DefaultSelector_FieldAzimuth, DefaultSelector_Latitude, &
    & DestroyL2PC, DestroyL2PCDatabase, DestroyBinSelectorDatabase, &
    & Diff, Dump, FlushL2PCBins, L2PC_T, &
    & LoadHessian, LoadMatrix, LoadVector, OutputHDF5L2PC, &
    & PopulateL2PCBin, PopulateL2PCBinByName, &
    & ReadCompleteHDF5L2PCFile

  interface DIFF
    module procedure DiffL2PCFiles, DiffL2PCs
  end interface

  interface DUMP
    module procedure DUMPONEL2PC, DUMPL2PCDATABASE, DUMPL2PCFILE
  end interface

  interface DUMP_PRIVATE
    module procedure DUMPONEL2PC, DUMPL2PCDATABASE, DUMPL2PCFILE, DUMPL2PCINFO
  end interface
  ! This is an update to the L2PCs where we can store both Jacobians and Hessians
  ! Previously the L2PC's were just Matrix_Ts, now they're more diversified
  ! As this contains pointers to vector_T's and so on, I maintain a private
  ! set of databases of these in this module.  We can't use the main databases,
  ! as these would get destroyed at the end of each chunk.

  ! The supporting databases
  integer, dimension(:), pointer, public, save :: FileIDDatabase => NULL()
  type(QuantityTemplate_T), dimension(:), pointer, save :: L2PCQTS => NULL()
  type(VectorTemplate_T), dimension(:), pointer, save :: L2PCVTS => NULL()
  type(Vector_T), dimension(:), pointer, save :: L2PCVS => NULL()

  logical :: verbose
  logical :: verboser

  ! This type holds an l2pc sometimes called a "bin"
  type L2PC_T
    integer :: NAME                     ! The name of the L2PC bin
    type(Matrix_T) :: J                 ! The Jacobian
    logical :: GOTH                     ! Set true if also have Hessian Info
    type(Hessian_T) :: H                ! The Hessian
  end type L2PC_T

  ! The L2PC database
  ! We could create much-needed flexibility if we allowed the user to maintain
  ! a separate database or multiple databases
  ! The same goes for binSelector and L2PCInfo databases
  type(L2PC_T), dimension(:), pointer, public, save :: L2PCDatabase => NULL()

  ! This datatype describes a selection rule for l2pc bins.
  type BinSelector_T
    integer :: selectorType             ! What quantity type does this apply to
    integer :: molecule                 ! What molecule does it apply to
    integer :: nameFragment             ! A possible name fragment
    logical :: exact                    ! We require an exact numeric match for this selector
    real(r8), dimension(2) :: heightRange ! The height range for this selector
    real(r8) :: cost                    ! The cost for that range
  end type BinSelector_T

  type(BinSelector_T), dimension(:), pointer, save :: BINSELECTORS => NULL()

  type L2PCInfo_T
    integer :: fileID     ! What is the HDF5 file ID
    integer :: binID      ! What is the groupID for the bin
    integer :: blocksID   ! What is the groupID for the blocks
    integer :: hBlocksID  ! What is the groupID for the hessian blocks
    integer, dimension(:,:), pointer :: BLOCKID => NULL()
    character(len=64) :: matrixName
  end type L2PCINFO_T
  type(L2PCInfo_T), dimension(:), pointer, save :: L2PCINFO => NULL()

  ! Default bin selectors (see CreateDefaultBinSelectors below)
  logical, parameter :: DEEBUG = .false.
  integer, parameter :: DEFAULTSELECTOR_LATITUDE = 1
  integer, parameter :: DEFAULTSELECTOR_FIELDAZIMUTH = 2
  
  logical, parameter :: DIEIFDESTROYFAILS = .false.
  integer, parameter :: MAXNBINS          = 200 ! max nBins in a file

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

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
    & call getString( STRING, STRING_TEXT, CAP, STRIP, NOERROR, IERR, &
    & START, END )
  end subroutine get_string

  ! ------------------------------------ AddBinSelectorToDatabase --
  integer function AddBinSelectorToDatabase ( database, item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    type (BinSelector_T), dimension(:), pointer :: DATABASE
    type (BinSelector_T) :: item
    ! Local variables
    type (BinSelector_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"
    AddBinSelectorToDatabase = newSize
  end function AddBinSelectorToDatabase

  ! ------------------------------------  Add fileID to database ----
  integer function AddFileIDToDatabase ( Database, Item )

    ! This function simply adds a fileID  to a database of this type

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    integer, dimension(:), pointer :: Database
    integer :: Item

    integer, dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddFileIDToDatabase = newSize
  end function AddFileIDToDatabase

  ! ------------------------------------  Add l2pc to database ----
  integer function AddL2PCToDatabase ( Database, Item )

    ! This function simply adds an l2pc  to a database of said l2pc s.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type(L2PC_T), dimension(:), pointer :: Database
    type(L2PC_T) :: Item

    type(L2PC_T), dimension(:), pointer :: TempDatabase

    include "addItemToDatabase.f9h"

    AddL2PCToDatabase = newSize
  end function AddL2PCToDatabase

  ! ----------------------------------- AdoptVectorTemplate --
  type (VectorTemplate_T) function AdoptVectorTemplate ( bin, name, quantityTemplates, source, message ) &
    & result ( vectorTemplate )
    ! This function is used to import a vector template from our private l2pc database
    ! into the mainstream database.  It also adopts the individual quantities into
    ! the mainstream quantity template database.

    ! Dummy arguments
    integer, intent(in) :: BIN          ! Name of MATRIX
    integer, intent(in) :: NAME         ! Name for new vector template
    type(QuantityTemplate_T), dimension(:), pointer :: QUANTITYTEMPLATES ! Mainstream database
    integer, intent(in) :: SOURCE       ! L_COLUMNS, L_ROWS
    character (len=*), intent(out) :: MESSAGE ! Error message

    ! Local variables
    integer :: I                        ! Array index
    type(VectorTemplate_T), pointer :: SOURCETEMPLATE ! Sought vector template
    type(Matrix_T), pointer :: M        ! The matrix
    integer :: QTY                      ! Loop counter
    type(QuantityTemplate_T) :: QT      ! A template

    ! Executable code
    message = ''
    call NullifyVectorTemplate ( vectorTemplate )

    i = FindFirst ( l2pcDatabase%name, bin )
    if ( i == 0 ) then
      message = 'No such l2pc bin found'
      return
    end if

    m => l2pcDatabase(i)%j
    if ( source == l_rows ) then
      sourceTemplate => m%row%vec%template
    else
      sourceTemplate => m%col%vec%template
    end if

    vectorTemplate = sourceTemplate
    vectorTemplate%name = name
    ! We need our own quantities array to point to the mainstream database
    nullify ( vectorTemplate%quantities )
    call Allocate_test ( vectorTemplate%quantities, vectorTemplate%noQuantities, &
      & 'adopted vectorTemplate%quantities', ModuleName )
    ! Copy each quantity from the l2pc database into the mainstream one
    do qty = 1, vectorTemplate%noQuantities
      call CopyQuantityTemplate ( qt, l2pcQTs ( sourceTemplate%quantities(qty) ) )
      i = FindFirst ( quantityTemplates%name, qt%name )
      if ( i /= 0 ) then
        ! Found this quantity already. Is it marked as one to be adopted
        if ( quantityTemplates(i)%quantityType == l_adopted ) then
          vectorTemplate%quantities(qty) = i
          ! Shallow copy this into the database over the 'adopted' placeholder,
          ! Destroy old placeholder first.
          call DestroyQuantityTemplateContents ( quantityTemplates(i) )
          quantityTemplates(i) = qt
        else
          call get_string ( qt%name, message, strip=.true. )
          message = 'Quantity ' // trim(message) // ' already exist, cannot adopt'
        end if
      else
        vectorTemplate%quantities(qty) = AddQuantityTemplateToDatabase ( quantityTemplates, qt )
      end if
      call NullifyQuantityTemplate ( qt ) ! To avoid treading on one in database
    end do
  end function AdoptVectorTemplate

  ! -------------------------------------- CreateDefaultBinSelectors --
  subroutine CreateDefaultBinSelectors
    ! This routine creates two bin selectors we know we need even
    ! if the user doesn't supply them in the l2cf.  If you add more
    ! here please make sure to also add DefaultSelector_... parmeters
    ! above.
    ! Local variables
    type (BinSelector_T) :: SEL
    integer :: DUMMY
    ! Executable code
    ! Create a latitude based one.  This is the one used by default if no
    ! other is supplied
    sel%selectorType = l_latitude
    sel%molecule = 0
    sel%nameFragment = 0
    sel%exact = .false.
    sel%heightRange = 0.0_r8
    sel%cost = 1.0_r8
    dummy = AddBinSelectorToDatabase ( binSelectors, sel )

    ! Create a template for a field_azimith one for the polarLinear model
    sel%selectorType = l_fieldAzimuth
    sel%molecule = 0
    sel%nameFragment = 0
    sel%exact = .true.
    sel%heightRange = 0.0_r8
    sel%cost = 1.0
    dummy = AddBinSelectorToDatabase ( binSelectors, sel )

  end subroutine CreateDefaultBinSelectors

  ! -------------------------------------- DestroyBinSelectorDatabase
  subroutine DestroyBinSelectorDatabase
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    ! Local variables
    integer(c_intptr_t) :: Addr            ! For tracing
    integer :: S, STATUS                   ! Flag from deallocate
    ! Executable code
    if ( .not. associated ( binSelectors ) ) return
    s = size(binSelectors) * storage_size(binSelectors) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(binSelectors(1)), addr)
    deallocate ( binSelectors, stat=status )
    call test_deallocate ( status, ModuleName, "binSelectors", s, address=addr )
  end subroutine DestroyBinSelectorDatabase

  ! ----------------------------------------------- DestroyL2PC ----
  subroutine DestroyL2PC ( l2pc )
    ! Dummy arguments
    type (L2pc_t), intent(inout), target :: L2PC

    integer :: QUANTITY                 ! Loop index
    integer :: VECTOR                   ! Loop index
    integer :: JH                       ! Loop index

    type (Vector_T), pointer :: V       ! Temporary pointer
    integer, parameter :: jacobians = 1
    integer, parameter :: hessians = jacobians + 1
    integer, parameter :: rows = 1
    integer, parameter :: cols = rows + 1

    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    do jh = jacobians, hessians
      if ( jh == hessians .and. .not. l2pc%gotH ) cycle
      if ( verbose ) &
        & call outputNamedValue( ' About to destroy 1 (j) or 2 (h)', jh )
      if ( jh == hessians ) then
        call output( 'Unable to destroy hessian quantity templates yet', advance='yes' )
        cycle
      endif
      do vector = rows, cols
        if ( verbose ) &
          & call outputNamedValue( ' About to destroy 1 (rows) or 2 (cols)', vector )
        if ( vector == cols ) then
          if ( jh == jacobians ) then
            v => l2pc%j%col%vec
          else
            v => l2pc%h%col%vec
          end if
        else
          if ( jh == jacobians ) then
            v => l2pc%j%row%vec
          else
            v => l2pc%h%row%vec
          end if
        end if
        if ( .not. associated(v) ) cycle
        if ( .not. associated(v%quantities) ) cycle
        if ( verbose ) then
          call output( ' About to destroy vectors', advance='yes' )
          call outputNamedValue( ' Num quantities', size(v%quantities) )
        endif
        do quantity = 1, size(v%quantities)
          call DestroyQuantityTemplateContents (v%quantities(quantity)%template )
        end do
        call DestroyVectorInfo ( v )
      end do
    end do

    ! Destroy kStar
    if ( verbose ) &
      & call output( ' About to destroy matrix', advance='yes' )
    call DestroyMatrix ( l2pc%j )
    if ( verbose ) &
      & call output( ' About to destroy hessian', advance='yes' )
    if ( l2pc%goth ) call DestroyHessian ( l2pc%h )

  end subroutine DestroyL2PC

  ! ------------------------------------------- DestroyL2PCDatabase ---
  subroutine DestroyL2PCDatabase

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, S, Status

    verbose = switchDetail( switches, 'l2pc') > -1
    if ( .not. associated(l2pcDatabase) ) then
      if ( verbose ) &
        & call output( 'l2pc database empty' )
    else
      do i = 1, size(l2pcDatabase)
        if ( verbose ) &
          & call outputNamedValue( 'Destroying l2pc db entry number ', i )
        call DestroyL2PC ( l2pcDatabase(i) )
      end do
      s = size(l2pcDatabase) * storage_size(l2pcDatabase) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(l2pcDatabase(1)), addr)
      deallocate ( l2pcDatabase, stat=status )
      call test_deallocate ( status, ModuleName, "l2pcDatabase", s, address=addr )
    end if

    ! Also destroy the info database (i.e. close files)
    call DestroyL2PCInfoDatabase
  end subroutine DestroyL2PCDatabase

  ! --------------------------------------- DiffL2PCFiles ---------------
  subroutine DiffL2PCFiles ( L2PCFile1, L2PCFile2, details, options )
    ! This subroutine Diffs an l2pc to stdout
    ! assuming it was read already
    use output_m, only: output

    ! Dummy arguments
    type (MLSFile_T), intent(inout) :: L2PCFile1, L2PCFile2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't Diff multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Diff even multi-dim arrays
    !                                        ! Default 0
    character(len=*), intent(in), optional :: options ! any of {jh*b[]}
                                             ! D: read files before Diffing
                                             ! j: Diff only matrices
                                             ! h: Diff only hessians
                                             ! *: Diff matrices and hessians
                                             ! b[HCN]: Diff only HCN blocks
                                             ! default is *
    ! Local variables
    logical :: different
    integer :: i
    integer, dimension(MAXNBINS) :: indices1
    integer, dimension(MAXNBINS) :: indices2
    logical :: mustRead
    integer :: nBins1
    integer :: nBins2
    logical :: silent
    logical :: verbose
    ! Executable
    mustRead = (optionDetail(options, 'D' ) == 'yes' )
    verbose = (optionDetail(options, 'v' ) == 'yes' )
    silent = (optionDetail(options, 'm' ) == 'yes' )
    if ( verbose ) call outputNamedValue( 'options', options )
    if ( verbose ) call outputNamedValue( 'mustRead', mustRead )
    if ( mustRead ) then
      call ReadCompleteHDF5L2PCFile ( L2PCFile1, 0, shallow=.true. )
      call ReadCompleteHDF5L2PCFile ( L2PCFile2, 0, shallow=.true. )
      if ( .not. associated(L2PCDataBase) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Sorry-l2pc database still empty' )
      endif
      do i=1, size(L2PCDataBase)
        if ( any( &
          & (/ L2PCFile1%FileID%f_id, L2PCFile2%FileID%f_id /) &
          & == fileIDDataBase(i) ) ) &
          & call PopulateL2PCBin( i )
      enddo
    endif
    call getNBinsInFile ( L2PCFile1, indices1, nBins1 )
    call getNBinsInFile ( L2PCFile2, indices2, nBins2 )
    if ( verbose ) call outputNamedValue( 'nBins1', nBins1 )
    if ( verbose ) call outputNamedValue( 'indices1', indices1 )
    if ( verbose ) call outputNamedValue( 'nBins2', nBins2 )
    if ( verbose ) call outputNamedValue( 'indices2', indices2 )
    different = .false.
    do i=1, min( nBins1, nBins2 )
      if ( .not. silent ) call newline
      if ( verbose ) call outputNamedValue( 'file index', i )
      call DiffL2PCs( L2PCDataBase(indices1(i)), L2PCDataBase(indices2(i)), &
        & details, options, different )
    enddo
    if ( different ) then
      call output( 'These L2PC files differ', advance='yes' )
    elseif ( silent ) then
      call output( 'These L2PC files do not differ', advance='yes' )
    endif
  end subroutine DiffL2PCFiles

  ! --------------------------------------- DiffL2PCs ---------------
  subroutine DiffL2PCs ( L2pc1, L2pc2, details, options, different )
    ! This subroutine diffs two l2pcs

    use hessianModule_1, only: Diff
    use matrixModule_1, only: Diff
    use MLSStrings, only: Delete

    ! Dummy arguments
    type (l2pc_t), intent(inout) :: L2pc1, l2pc2
    integer, intent(in), optional :: DETAILS ! passed to Vector and Matrix dumps
    character(len=*), intent(in), optional :: options ! any of {jh*b[]}
                                             ! j: dump only matrices
                                             ! h: dump only hessians
                                             ! *: dump matrices and hessians
                                             ! b[HCN]: dump only HCN blocks
                                             ! default is *
    logical, intent(inout), optional :: different ! Are matrices different?
    ! Local variables
    character(len=16) :: audible
    logical :: dumpHessians
    logical :: dumpJacobians
    logical :: itsDifferent
    integer :: myDetails
    logical :: silent                        ! if opttions contain 'm'
    logical :: verbose                       ! if opttions contain 'b'
    ! Executable code
    myDetails = 0
    if ( present(details) ) myDetails = details
    dumpHessians = .true.
    dumpJacobians = .true.
    verbose = .false.
    silent = .false.
    audible = ' '
    if ( present(options) ) then
      dumpHessians = optionDetail( options, 'h' ) == 'yes'
      dumpJacobians = optionDetail( options, 'j' ) == 'yes'
      verbose = (optionDetail(options, 'v' ) == 'yes' )
      silent = (optionDetail(options, 'm' ) == 'yes' )
      audible = Delete ( options, 'm' )
      if ( verbose ) then
        call outputNamedValue( 'options', trim(options) )
        call outputNamedValue( 'dumpJacobians', dumpJacobians )
        call outputNamedValue( 'dumpHessians', dumpHessians )
      endif
    endif
    if ( verbose ) call outputNamedValue( 'myDetails', myDetails )
    if ( .not. silent ) then
      call output( '- Diffing L2PCs -', advance='yes' )
      call NameThoseBins
      if ( verbose ) call outputNamedValue( 'l2pc1%goth', l2pc1%goth )
      if ( verbose ) call outputNamedValue( 'l2pc2%goth', l2pc2%goth )
    endif
    itsDifferent = .false.
    if ( dumpJacobians ) then
      call diff ( l2pc1%j, l2pc2%j, &
        & details, options=options, different=itsDifferent )
      if ( silent .and. itsDifferent ) &
        & call diff ( l2pc1%j, l2pc2%j, details, options=audible )
    endif

    if ( l2pc1%goth .and. dumpHessians ) &
      & call diff ( l2pc1%h, l2pc2%h, details, options=options )
    if ( present(different) ) different = different .or. itsDifferent
    if ( silent .and. .not. itsDifferent ) then
      call NameThoseBins
      call output( 'These L2PC bins do not differ', advance='yes' )
    endif
  contains
    subroutine NameThoseBins
      if ( l2pc1%name > 0 ) then
        call display_string ( l2pc1%name, before='name: ' )
        call output ( l2pc1%name, before=' (', after=')', advance='yes' )
      else
        call output( '*** Uh-oh, name not found in string table', advance='yes' )
      endif
      if ( l2pc2%name > 0 ) then
        call display_string ( l2pc2%name, before='name: ' )
        call output ( l2pc2%name, before=' (', after=')', advance='yes' )
      else
        call output( '*** Uh-oh, name not found in string table', advance='yes' )
      endif

    end subroutine NameThoseBins
  end subroutine DiffL2PCs

  ! --------------------------------------- DumpL2PCDatabase ---------------
  subroutine DumpL2PCDatabase ( L2pcDB, details, onlytheseblocks, options )
    ! This subroutine dumps an l2pc to stdout

    ! Dummy arguments
    type (l2pc_t), dimension(:), intent(in), target :: L2pcDB
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    character(len=*), dimension(3), intent(in), optional :: ONLYTHESEBLOCKS ! passed to Vector and Matrix dumps
    character(len=*), intent(in), optional :: options ! any of {mh*b[]}
                                             ! m: dump only matrices
                                             ! h: dump only hessians
                                             ! v: dump only x*, y* vectors
                                             ! *: dump everything
                                             ! b[HCN]: dump only HCN blocks
                                             ! default is *
    ! Local variables
    integer :: i
    ! Executable
    do i=1, size(L2PCDB)
      call DumpOneL2PC( L2PCDB(i), details, ONLYTHESEBLOCKS, options )
    enddo
  end subroutine DumpL2PCDatabase

  ! --------------------------------------- DumpL2PCFile ---------------
  subroutine DumpL2PCFile ( L2PCFile, details, onlytheseblocks, options )
    ! This subroutine dumps an l2pc to stdout
    ! assuming it was read already

    ! Dummy arguments
    type (MLSFile_T), intent(inout) :: L2PCFile
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    character(len=*), dimension(3), intent(in), optional :: ONLYTHESEBLOCKS ! passed to Vector and Matrix dumps
    character(len=*), intent(in), optional :: options ! any of {mh*b[]}
                                             ! r: read file before dumping
                                             ! m: dump only matrices
                                             ! h: dump only hessians
                                             ! v: dump only x*, y* vectors
                                             ! *: dump everything
                                             ! b[HCN]: dump only HCN blocks
                                             ! default is *
    ! Local variables
    integer :: i
    logical :: mustRead
    ! Executable
    mustRead = .false.
    if ( present(options) ) mustRead = (optionDetail(options, 'r' ) == 'yes' )
    if ( mustRead ) then
      call ReadCompleteHDF5L2PCFile ( L2PCFile, 0, shallow=.true. )
      do i=1, size(L2PCDataBase)
        if ( L2PCFile%FileID%f_id == fileIDDataBase(i) ) &
          & call PopulateL2PCBin( i )
      enddo
    endif
    call dumpMLSFile( L2PCFile, details=1 )
    call dump ( fileIDDataBase, 'file id database', format='(i10)' )
    do i=1, size(L2PCDataBase)
      if ( L2PCFile%FileID%f_id /= fileIDDataBase(i) ) cycle
      call newline
      call outputNamedValue( 'db index', i )
      call DumpOneL2PC( L2PCDataBase(i), details, onlytheseblocks, options )
    enddo
  end subroutine DumpL2PCFile

  ! --------------------------------------- DumpOneL2PC ---------------
  subroutine DumpOneL2PC ( L2pc, details, ONLYTHESEBLOCKS, options )
    ! This subroutine dumps an l2pc to stdout

    use hessianModule_1, only: dump, dump_layout
    use MatrixModule_1, only: dump_layout

    ! Dummy arguments
    type (l2pc_t), intent(in), target :: L2pc
    integer, intent(in), optional :: DETAILS ! passed to Vector and Matrix dumps
    character(len=*), dimension(3), intent(in), optional :: ONLYTHESEBLOCKS ! passed to Vector and Matrix dumps
    character(len=*), intent(in), optional :: options ! any of {Lmh*b[]d[]}
                                             ! L: Just show the top level layouts
                                             ! m: dump only matrices
                                             ! h: dump only hessians
                                             ! v: dump only x*, y* vectors
                                             ! *: dump everything
                                             ! b[HCN]: dump only HCN blocks
                                             ! d[dopts]: pass dopts when dumping arrays
                                             ! default is *
    ! Local variables
    logical :: dumpHessians
    logical :: dumpJacobians
    logical :: dumpVectors
    integer :: i
    logical :: LayoutOnly
    integer :: myDetails
    ! character(len=128) :: myMolecules

    ! Executable code
    myDetails = 0
    if ( present(details) ) myDetails = details
    dumpHessians = .true.
    dumpJacobians = .true.
    dumpVectors = .true.
    if ( present(options) ) then
      dumpHessians = optionDetail( options, 'h' ) == 'yes'
      dumpJacobians = optionDetail( options, 'm' ) == 'yes'
      dumpVectors = optionDetail( options, 'v' ) == 'yes'
      LayoutOnly = optionDetail( options, 'L' ) == 'yes'
    endif
    ! call outputNamedValue( 'myDetails', myDetails )
    call output( '- Dump of L2PC -', advance='yes' )
    if ( l2pc%name > 0 ) then
      call display_string ( l2pc%name, before='name: ' )
      call output ( l2pc%name, before=' (', after=')', advance='yes' )
    else
      call output( "(name not found in string table)", advance='yes' )
    endif
    ! First dump the xStar and yStar
    if ( dumpVectors ) then
      if ( LayoutOnly ) then
        call output ('Layout(xStar)', advance='yes' )
        do i=1, size(l2pc%j%col%vec%quantities)
          call display_string ( l2pc%j%col%vec%quantities(i)%template%name, &
          & strip=.true., advance='yes' )
        enddo
        call output ('Layout(yStar)', advance='yes' )
        do i=1, size(l2pc%j%row%vec%quantities)
          call display_string ( l2pc%j%row%vec%quantities(i)%template%name, &
          & strip=.true., advance='yes' )
        enddo
      else
        call dump ( l2pc%j%col%vec, details=details, name='xStar' )
        call dump ( l2pc%j%row%vec, details=details, name='yStar' )
      endif
    endif

    if ( dumpJacobians ) then
      ! Now dump kStar
      if ( LayoutOnly ) then
        call dump_layout( l2pc%j, 'kStar' )
      else
        call dump ( l2pc%j, 'kStar', details )
      endif
    endif

    ! Now dump the Hessian
    if ( l2pc%goth .and. dumpHessians ) then
      if ( layoutOnly ) then
        call dump_layout( l2pc%h, 'hStar' )
      else
        call dump ( l2pc%h, 'hStar', details, onlyTheseBlocks, options )
      endif
    endif
  end subroutine DumpOneL2PC

  ! --------------------------------------- DumpL2PCInfo ---------------
  subroutine DumpL2PCInfo ( L2PCInfo, details )
    ! This subroutine dumps an l2pc to stdout

    ! Dummy arguments
    type (L2PCInfo_T), intent(in) :: L2PCInfo
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but size
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    ! Local variables
    integer :: myDetails
    ! Executable
    myDetails = 0
    if ( present(details) ) myDetails = details
    call outputNamedValue ( 'fileID', L2PCInfo%fileID )
    call outputNamedValue ( 'binID', L2PCInfo%binID )
    call outputNamedValue ( 'blocksID', L2PCInfo%blocksID )
    call outputNamedValue ( 'hblocksID', L2PCInfo%hBlocksID )
    call outputNamedValue ( 'matrixName', trim(L2PCInfo%matrixName) )
    if ( myDetails > 0 ) call dump( L2PCInfo%blockID, 'BlockID' )
  end subroutine DumpL2PCInfo

  ! ------------------------------------ FlushL2PCBins -------------
  subroutine FlushL2PCBins
    ! Local variables
    integer :: BIN              ! Loop counter
    integer :: BLOCKROW         ! Loop counter
    integer :: BLOCKCOL         ! Loop counter
    integer :: I, J, K

    type ( L2PC_T ), pointer          :: L2PC
    type ( MatrixElement_T), pointer  :: M0
    type ( HessianElement_T), pointer :: H0

    ! Executable code
    if ( .not. associated ( l2pcDatabase ) ) return
    do bin = 1, size ( l2pcDatabase )
      l2pc => l2pcDatabase ( bin )
      do blockRow = 1, l2pc%j%row%NB
        do blockCol = 1, l2pc%j%col%NB
          m0 => l2pc%j%block ( blockRow, blockCol )
          if ( m0%kind /= m_absent ) then
            call DestroyBlock ( m0 )
            m0%kind = M_Unknown
          end if
        end do
      end do
    end do
    if ( l2pc%goth ) then
      do i = 1, l2pc%h%row%NB
        do j = 1, l2pc%h%col%NB
          do k = 1, l2pc%h%col%NB
            h0 => l2pc%h%block(i,j,k)
            if ( h0%kind /= h_absent ) then
              call DestroyBlock ( h0 )
              h0%kind = h_unknown
            end if
          end do
        end do
      end do
    endif
  end subroutine FlushL2PCBins

  ! ---------------------------------- LoadHessian ----------
  subroutine LoadHessian ( Hessian, name, message )
    ! This function copies the Hessian from the l2pc database into the
    ! given Hessian (presumably from the main dabase)

    ! Dummy arguments
    type(Hessian_T), intent(inout) :: Hessian
    integer, intent(in) :: NAME
    character (len=*), intent(out) :: MESSAGE

    ! Local variables
    integer :: I                        ! Index
    type(Hessian_T), pointer :: SOURCE

    ! Executable code
    message = ''
    i = FindFirst ( l2pcDatabase%name, name )
    call outputNamedValue( 'l2pcdatabase names', l2pcDatabase%name )
    call outputNamedValue( 'bin name', name )
    call outputNamedValue( 'i', i )
    if ( i == 0 ) then
      message = 'No such l2pc Hessian with your name'
      call display_string( name, advance='yes' )
      return
    end if
    if ( .not. l2pcDatabase(i)%goth ) then
      message = 'This bin has no Hessian to copy'
      return
    endif
    call PopulateL2PCBin ( i, ignoreHessian=.false. )
    source => l2pcDatabase(i)%h
    if ( .not. DoVectorsMatch ( Hessian%row%vec, source%row%vec, verbose=.true. ) ) then
      message = 'Rows do not match for loading'
      call output( 'Hessian row vector', advance='yes' )
      call dump( Hessian%row%vec, details=0 )
      call output( 'source row vector', advance='yes' )
      call dump( source%row%vec, details=0 )
      return
    endif
    if ( .not. DoVectorsMatch ( Hessian%col%vec, source%col%vec ) ) then
      message = 'Columns do not match for loading'
      return
    endif
    call CopyHessianValue ( Hessian, source, allowNameMismatch=.true. )
  end subroutine LoadHessian

  ! ---------------------------------- LoadMatrix ----------
  subroutine LoadMatrix ( matrix, name, message, IgnoreHessian )
    ! This function copies from source in the l2pc database into a
    ! specified matrix
    ! It checks that they are compatible according to many criteria (too many?)

    ! Dummy arguments
    type(Matrix_T), intent(inout) :: MATRIX
    integer, intent(in) :: NAME
    character (len=*), intent(out) :: MESSAGE
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    integer :: I                        ! Index
    type(Matrix_T), pointer :: SOURCE

    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    message = ''
    i = FindFirst ( l2pcDatabase%name, name )
    call outputNamedValue( 'l2pcdatabase names', l2pcDatabase%name )
    call outputNamedValue( 'bin name', name )
    call outputNamedValue( 'i', i )
    if ( i == 0 ) then
      message = 'No such l2pc matrix'
      return
    end if
    call PopulateL2PCBin ( i, ignoreHessian )
    source => l2pcDatabase(i)%j
    if ( .not. DoVectorsMatch ( matrix%row%vec, source%row%vec, verbose ) ) then
      message = 'Rows do not match for loading'
      call output( 'matrix row vector', advance='yes' )
      call dump( matrix%row%vec, details=0 )
      call output( 'source row vector', advance='yes' )
      call dump( source%row%vec, details=0 )
      return
    endif
    if ( .not. DoVectorsMatch ( matrix%col%vec, source%col%vec, verbose ) ) then
      message = 'Columns do not match for loading'
      call output( 'matrix col vector', advance='yes' )
      call dump( matrix%col%vec, details=0 )
      call output( 'source col vector', advance='yes' )
      call dump( source%col%vec, details=0 )
      return
    endif
    call CopyMatrixValue ( matrix, source, allowNameMismatch=.true. )
  end subroutine LoadMatrix

  ! ---------------------------------- LoadVector ----------
  subroutine LoadVector ( vector, name, source, message, IgnoreHessian )
    ! This function copies the matrix from the l2pc database into the
    ! given matrix (presumably from the main dabase)

    ! Dummy arguments
    type(Vector_T), intent(inout) :: VECTOR
    integer, intent(in) :: NAME
    integer, intent(in) :: SOURCE       ! l_rows or l_columns
    character (len=*), intent(out) :: MESSAGE
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    integer :: I                        ! Index
    type(Matrix_T), pointer :: MATRIX
    type(Vector_T), pointer :: SOURCEVECTOR

    ! Executable code
    message = ''
    i = FindFirst ( l2pcDatabase%name, name )
    if ( i == 0 ) then
      message = 'No such l2pc matrix'
      return
    end if
    call PopulateL2PCBin ( i, ignoreHessian )
    matrix => l2pcDatabase(i)%j
    select case ( source )
    case ( l_rows )
      sourceVector => matrix%row%vec
    case ( l_columns )
      sourceVector => matrix%col%vec
    end select
    if ( .not. DoVectorsMatch ( vector, sourceVector ) ) then
      message = 'Vector does not match for loading'
      return
    endif
    ! Use the clone option to avoid fussyness with vector templates
    call CopyVector ( vector, sourceVector, allowNameMismatch=.true. )
  end subroutine LoadVector

  ! --------------------------------------------- OutputHDF5L2PC
  subroutine OutputHDF5L2PC ( filename, matrices, hessians, &
    & quantitiesNode, secondDerivNode, packed, dontPack )
  use hdf5, only: h5fcreate_f, h5fclose_f, h5f_acc_trunc_f
  use MLSSTringLists, only: catLists
  use tree, only: sub_rosa
    character (len=*), intent(in) :: FILENAME
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES
    type (Hessian_T), dimension(:), pointer :: HESSIANS
    integer, intent(in) :: QUANTITIESNODE
    integer, intent(in) :: SECONDDERIVNODE
    logical, intent(in) :: PACKED
    integer, dimension(:), pointer :: DONTPACK

    ! Local variables
    integer :: DB_INDEX                 ! Index of matrix
    integer :: FIELD                    ! Node index
    integer :: FILEID                   ! ID of file
    logical :: GOTH                     ! True if there is a Hessian
    logical :: GOTM                     ! True if there is a Jacobian
    integer :: J
    integer :: NXT_INDEX                ! Index of next matrix (Hessian?)
    integer :: STATUS                   ! From HDF
    integer :: THISMOLECULE
    character(len=32) :: moleculeName
    character(len=32), dimension(3) :: ONLYTHESEBLOCKS ! passed to Vector and Matrix dumps
    type (Matrix_T), pointer :: tmpMatrix
    type (Hessian_T), pointer :: tmpHessian

    ! Executable code
    call H5FCreate_F ( trim(filename), H5F_ACC_TRUNC_F, fileID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 l2pc file for output.' )
    ONLYTHESEBLOCKS = '*'
    if ( SECONDDERIVNODE > 0 ) then
      ONLYTHESEBLOCKS(2) = ' '
      do j = 2, nsons(secondDerivNode)
        thisMolecule = decoration( subtree( j, secondDerivNode ) )
        call get_string( sub_rosa(thisMolecule), moleculeName )
        onlyTheseBlocks(2) = catLists( onlyTheseBlocks(2), thisMolecule )
      end do
      onlyTheseBlocks(3) = onlyTheseBlocks(2)
    endif
    nullify ( tmpMatrix, tmpHessian )
    field = 2
    do
      if ( field > nsons(quantitiesNode) ) exit
      db_index = decoration(decoration(subtree(field, quantitiesNode )))
      if ( associated(matrices) ) then
        call GetActualMatrixFromDatabase ( matrices(db_index), tmpMatrix )
        gotm = .true.
      endif
      goth = .false.
      if ( field < nsons(quantitiesNode) ) then
        nxt_index = decoration(decoration(subtree(field+1, quantitiesNode )))
        if ( nxt_index < 0 ) then
          tmpHessian => hessians ( -nxt_index )
          goth = .true.
          field = field + 1
        end if
      end if
      if ( goth ) then
        call writeOneHDF5L2PC ( tmpMatrix, fileID, packed, dontPack, &
          & hessian=tmpHessian, onlyTheseBlocks=onlyTheseBlocks )
      elseif ( gotm ) then
        call writeOneHDF5L2PC ( tmpMatrix, fileID, packed, dontPack )
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to write any Jacobians or Hessians to this l2pc file.' )
      end if
      field = field + 1
    end do ! Loop over fields
    call H5FClose_F ( fileID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Unable to close hdf5 l2pc file.' )

  end subroutine OutputHDF5L2PC

  ! --------------------------------------- WriteOneHDF5L2PC -----------
  subroutine WriteOneHDF5L2PC ( JACOBIAN, fileID, packed, dontPack, &
    & hessian, onlyTheseBlocks )
    use hessianModule_0, only: optimizeBlock
    use HDF5, only: h5gclose_f, h5gcreate_f
    use MLSKinds, only: rm
    use MLSHdf5, only: makeHDF5Attribute, saveAsHDF5DS
    ! This subroutine writes an l2pc to a file in hdf5 format

    ! Dummy arguments
    type (matrix_T), intent(in), target :: JACOBIAN
    integer, intent(in) :: fileID
    logical, intent(in) :: PACKED
    integer, dimension(:), pointer :: DONTPACK
    type (Hessian_T), intent(in), optional :: HESSIAN
    character(len=*), dimension(3), intent(in), optional :: ONLYTHESEBLOCKS ! passed to Vector and Matrix dumps

    ! Local variables
    integer :: I, J, K                  ! Loop index (row, col, col)
    integer :: BlocksGID                ! ID of group containing all blocks
    integer :: BlockGID                 ! ID of this block group
    integer :: matrixID                 ! ID of hdf5 group containing matrix
    integer :: STATUS                   ! Error flag

    integer, dimension(jacobian%row%nb) :: ROWBLOCKMAP
    integer, dimension(jacobian%col%nb) :: COLBLOCKMAP
    logical, dimension(jacobian%row%vec%template%noQuantities), target :: ROWPACK
    logical, dimension(jacobian%col%vec%template%noQuantities), target :: COLPACK

    character ( len=32 ) :: NAME        ! A name for output

    type (MatrixElement_T), pointer :: M0 ! A Matrix_0 within kStar
    type (HessianElement_T), pointer :: H0 ! A Hessian_0 within hStar
    logical :: verbose

    ! Executable code
    verbose = ( switchDetail(switches, 'hess') > -1 ) .or. DEEBUG
    ! If we have a Hessian, make sure it matches the jacobian
    if ( present ( hessian ) ) then
      if ( jacobian%row%vec%name /= hessian%row%vec%name ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Rows of Jacobian and Hessian do not match' )
      if ( jacobian%col%vec%name /= hessian%col%vec%name ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Columns of Jacobian and Hessian do not match' )
      if ( jacobian%row%instFirst .neqv. hessian%row%instFirst ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Row order of Jacobian and Hessian do not match' )
      if ( jacobian%col%instFirst .neqv. hessian%col%instFirst ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Column order of Jacobian and Hessian do not match' )
      if ( verbose ) then
        call get_string ( hessian%name, name, strip=.true., cap=.true. )
        call output( 'Writing Hessian named: ' // trim(name), advance='yes' )
        call outputNamedValue( 'Hesssian optimized already', hessian%optimizedAlready )
      endif
    end if

    ! Work out which quantities we can skip
    if ( packed ) then
      call MakeMatrixPackMap ( jacobian, rowPack, colPack, rowBlockMap, colBlockMap, dontPack, hessian )
    else
      rowPack = .true.
      colPack = .true.
      do i = 1, jacobian%row%nb
        rowBlockMap(i) = i
      end do
      do j = 1, jacobian%col%nb
        colBlockMap(j) = j
      end do
    end if

    ! First create the group for this.
    call get_string ( jacobian%name, name, strip=.true., cap=.true. )
    call h5gCreate_f ( fileID, trim(name), matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for l2pc matrix' )

    ! Now dump the vectors
    call WriteVectorAsHDF5L2PC ( matrixID, jacobian%col%vec, 'Columns', colPack )
    call WriteVectorAsHDF5L2PC ( matrixID, jacobian%row%vec, 'Rows', rowPack )

    ! Now create a group for the blocks
    call h5gCreate_f ( matrixID, 'Blocks', blocksGID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for all l2pc matrix blocks' )

    ! Make the flags
    call MakeHDF5Attribute ( blocksGID, 'rowInstanceFirst', jacobian%row%instFirst )
    call MakeHDF5Attribute ( blocksGID, 'colInstanceFirst', jacobian%col%instFirst )

    ! Now loop over the blocks and write them.
    do i = 1, jacobian%row%NB
      do j = 1, jacobian%col%NB
        ! Do we write this block?
        if ( rowPack(jacobian%row%quant(i)) .and. &
          &  colPack(jacobian%col%quant(j)) ) then
          ! Identify the block
          m0 => jacobian%block(i, j)
          ! Get a name for this group for the block
          ! write ( name, * ) 'Block', rowBlockMap(i), colBlockMap(j)
          call CreateBlockName ( rowBlockMap(i), colBlockMap(j), 0, name )
          ! Create a grop for this block
          call h5gCreate_f ( blocksGID, trim(name), blockGID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to create group for jacobian matrix block' )
          ! Stick some attributes on it
          call get_string ( &
            & jacobian%row%vec%quantities(&
            &    jacobian%row%quant(i))%template%name, name )
          call MakeHDF5Attribute ( blockGID, 'rowQuantity', trim(name) )
          call get_string ( &
            & jacobian%col%vec%quantities(&
            &    jacobian%col%quant(j))%template%name, name )
          call MakeHDF5Attribute ( blockGID, 'colQuantity', trim(name) )
          call MakeHDF5Attribute ( blockGID, 'rowInstance', jacobian%row%inst(i) )
          call MakeHDF5Attribute ( blockGID, 'colInstance', jacobian%col%inst(j) )
          call MakeHDF5Attribute ( blockGID, 'kind', m0%kind )
          ! Write the datasets
          if ( m0%kind /= m_absent ) then
            call SaveAsHDF5DS ( blockGID, 'values', real ( m0%values, r4 ) )
            if ( m0%kind /= m_full ) then
              call MakeHDF5Attribute ( blockGID, 'noValues', size(m0%values) )
              call SaveAsHDF5DS ( blockGID, 'r1', m0%r1 )
              call SaveAsHDF5DS ( blockGID, 'r2', m0%r2 )
            end if                      ! Not sparse or banded
          end if                        ! Not absent
          ! Close group for block
          call h5gClose_f ( blockGID, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to close group for l2pc matrix block' )
        end if                          ! Do this block
      end do
    end do

    ! Now close blocks group
    call h5gClose_f ( blocksGID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close qroup for all l2pc matrix blocks' )

    ! Create a group for the Hessian blocks if any
    if ( present ( hessian ) ) then
      call h5gCreate_f ( matrixID, 'HessianBlocks', blocksGID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create group for all l2pc matrix hessian blocks' )

      ! Make the flags
      call MakeHDF5Attribute ( blocksGID, 'rowInstanceFirst', hessian%row%instFirst )
      call MakeHDF5Attribute ( blocksGID, 'colInstanceFirst', hessian%col%instFirst )

      ! Now loop over the blocks and write them.
      do i = 1, hessian%row%NB
        do j = 1, hessian%col%NB
          do k = 1, hessian%col%NB
            ! Do we write this block?
            if ( rowPack(hessian%row%quant(i)) .and. &
              &  colPack(hessian%col%quant(j)) .and. &
              &  colPack(hessian%col%quant(k)) ) then
              ! Identify the block
              h0 => hessian%block ( i, j, k )
              if ( verbose ) then
                call outputNamedValue( 'block optimized already?', h0%optimizedAlready )
              endif
              call optimizeBlock ( h0 )
              ! Change kind of Hessian to absent if number of nonzero values is 0
              if ( h0%kind == h_sparse .and. &
                & (h0%tuplesFilled < 1 .or. .not. associated(h0%tuples) ) &
                & ) &
                & h0%kind = h_absent
              if ( h0%kind == h_full ) then
                if ( associated(h0%values) ) then
                  if ( count(h0%values /= 0._rm) < 1 ) h0%kind = h_absent
                else
                  h0%kind = h_absent
                endif
              endif
              ! Get a name for this group for the block
              call CreateBlockName ( rowBlockMap(i), colBlockMap(j), colBlockMap(k), name )
              if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call outputNamedValue( 'Block name', trim(name) )
              ! Create a group for this block
              call h5gCreate_f ( blocksGID, trim(name), blockGID, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to create group for hessian matrix block' )
              ! Stick some attributes on it
              ! Row name
              call get_string ( &
                & hessian%row%vec%quantities(&
                &    hessian%row%quant(i))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'rowQuantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'rowInstance', hessian%row%inst(i) )
              if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call outputNamedValue( 'Row name', trim(name) )
              ! First column name
              call get_string ( &
                & hessian%col%vec%quantities(&
                &    hessian%col%quant(j))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'col1Quantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'col1Instance', hessian%col%inst(j) )
              if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call outputNamedValue( '1st column name', trim(name) )
              ! Second column name
              call get_string ( &
                & hessian%col%vec%quantities(&
                &    hessian%col%quant(k))%template%name, name )
              call MakeHDF5Attribute ( blockGID, 'col2Quantity', trim(name) )
              call MakeHDF5Attribute ( blockGID, 'col2Instance', hessian%col%inst(j) )

              call MakeHDF5Attribute ( blockGID, 'kind', h0%kind )
              if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call outputNamedValue( '2nd column name', trim(name) )
              ! Write the datasets
              select case ( h0%kind )
              case ( h_full )
                call SaveAsHDF5DS ( blockGID, 'values', real ( h0%values, r4 ) )
              case ( h_sparse )
                call MakeHDF5Attribute ( blockGID, 'noValues', h0%tuplesFilled )
                call SaveAsHDF5DS ( blockGID, 'i', h0%tuples(1:h0%tuplesFilled)%i )
                call SaveAsHDF5DS ( blockGID, 'j', h0%tuples(1:h0%tuplesFilled)%j )
                call SaveAsHDF5DS ( blockGID, 'k', h0%tuples(1:h0%tuplesFilled)%k )
                call SaveAsHDF5DS ( blockGID, 'h', h0%tuples(1:h0%tuplesFilled)%h )
                if ( switchDetail( switches, 'hess' ) > 0 ) then
                  call outputNamedValue( 'noValues', h0%tuplesFilled )
                  call dump( h0%tuples(1:h0%tuplesFilled)%h, 'h)%tuple values' )
                end if
              end select
              ! Close group for block
              call h5gClose_f ( blockGID, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to close group for l2pc matrix block' )
            end if                          ! Do this block
          end do
        end do
      end do

      ! Now close blocks group
      call h5gClose_f ( blocksGID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close qroup for all l2pc matrix blocks' )
    end if

    ! Close matrix group
    call h5gClose_f ( matrixID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group for l2pc matrix' ) 
 end subroutine WriteOneHDF5L2PC

  ! ======================================= PRIVATE PROCEDURES ====================

  ! ------------------------------------ AddL2PCInfoToDatabase --
  integer function AddL2PCInfoToDatabase ( database, item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    type (L2PCInfo_T), dimension(:), pointer :: DATABASE
    type (L2PCInfo_T) :: item
    ! Local variables
    type (L2PCInfo_T), dimension(:), pointer :: TEMPDATABASE

    include "addItemToDatabase.f9h"
    AddL2PCInfoToDatabase = newSize
  end function AddL2PCInfoToDatabase

  ! ------------------------------------ DestroyL2PCInfoDatabase ----
  subroutine DestroyL2PCInfoDatabase

    use Allocate_Deallocate, only: Test_Deallocate
    use HDF5, only: H5FClose_f, H5GClose_f, H5ESet_Auto_f
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Local variables
    integer(c_intptr_t) :: Addr ! For tracing
    integer :: I                ! Loop counter
    integer :: S                ! Size in bytes of an object to deallocate
    integer :: STATUS           ! Flag from HDF

    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    if ( .not. associated(l2pcInfo) ) then
      if ( verbose ) &
        & call output( 'l2pcinfo database empty' )
      return
    endif
    call h5eSet_auto_f ( 0, status )
    do i = 1, size ( l2pcInfo )
      if ( verbose ) &
        & call outputNamedValue( 'Destroying l2pcinfo db entry number ', i )
      call h5gClose_f ( l2pcInfo(i)%blocksID, status )
      if ( status /= 0 .and. DIEIFDESTROYFAILS ) then
        call dump_private( l2pcInfo(i) )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close Blocks group for preserved input l2pc' )
      endif
      if ( l2pcInfo(i)%hblocksID /= 0 ) then
        call h5gClose_f ( l2pcInfo(i)%hblocksID, status )
        if ( status /= 0 .and. DIEIFDESTROYFAILS ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close Hessian Blocks group for preserved input l2pc' )
      end if
      call h5gClose_f ( l2pcInfo(i)%binID, status )
      if ( status /= 0 .and. DIEIFDESTROYFAILS ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close matrix group for preserved input l2pc' )
      ! Close the file?
      if ( count ( l2pcInfo%fileID == l2pcInfo(i)%fileID ) == 1 .and. &
        & all(l2pcInfo%fileID > 0 ) ) then
        ! We're the only one (left?) with this file, close it.
        if ( verbose ) &
          & call outputNamedValue('Closing hdf5 fileID', l2pcInfo(i)%fileID )
        call h5fClose_f ( l2pcInfo(i)%fileID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close hdf5 preserved input l2pc file' )
      end if
      l2pcInfo(i)%fileID = 0
    end do
    call h5eSet_auto_f ( 1, status )
    s = size(l2pcInfo) * storage_size(l2pcInfo) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(l2pcInfo(1)), addr)
    deallocate ( l2pcInfo, stat=i )
    call test_deallocate ( i, ModuleName, 'l2pcInfo', s, address=addr )

  end subroutine DestroyL2PCInfoDatabase

  ! ----------------------------------- MakeMatrixPackMap -----------
  subroutine MakeMatrixPackMap ( j, rowPack, colPack, rowBlockMap, colBlockMap, dontpack, h )
    ! This subroutine fills the boolean arrays rowPack, colPack
    ! (each length row/col%noQuantities) with a flag set true
    ! if the quantity has any derivatives at all

    ! Dummy arguments
    type (Matrix_t), intent(in) :: J
    logical, intent(out), dimension(J%row%vec%template%noQuantities) :: ROWPACK
    logical, intent(out), dimension(J%col%vec%template%noQuantities) :: COLPACK
    integer, intent(out), dimension(J%row%NB) :: ROWBLOCKMAP
    integer, intent(out), dimension(J%col%NB) :: COLBLOCKMAP
    integer, dimension(:), pointer :: DONTPACK ! Quantities not to pack
    type (Hessian_T), intent(in), optional :: H ! A hessian to consider also

    ! Local variables
    integer :: ROWQ                     ! Loop counter
    integer :: COLQ                     ! Loop counter
    integer :: COLQ1                    ! Loop counter
    integer :: COLQ2                    ! Loop counter
    integer :: ROWI                     ! Loop counter
    integer :: COLI                     ! Loop counter
    integer :: COLI1                    ! Loop counter
    integer :: COLI2                    ! Loop counter
    integer :: ROWBLOCK                 ! Block index
    integer :: COLBLOCK                 ! Block index
    integer :: COLBLOCK1                ! Block index
    integer :: COLBLOCK2                ! Block index
    logical, dimension(J%row%NB) :: ROWBLOCKFLAG ! Flags per row block
    logical, dimension(J%col%NB) :: COLBLOCKFLAG ! Flags per row block

    ! Executable code
    rowPack = .false.
    colPack = .false.

    ! Do a nested loop over cols/rows
    ! I tried to be fancy with cycles etc. but the code got really messy.
    ! This simple approach is probably the clearest.
    do rowQ = 1, j%row%vec%template%noQuantities
      do rowI = 1, j%row%vec%quantities(rowQ)%template%noInstances
        rowBlock = FindBlock ( j%row, rowQ, rowI )
        do colQ = 1, j%col%vec%template%noQuantities
          do colI = 1, j%col%vec%quantities(colQ)%template%noInstances
            colBlock = FindBlock ( j%col, colQ, colI )
            if ( j%block ( rowBlock, colBlock ) % kind /= M_Absent ) then
              rowPack ( rowQ ) = .true.
              colPack ( colQ ) = .true.
            end if
          end do
        end do
      end do
    end do

    ! Do the same thing for any hessian supplied
    if ( present ( h ) ) then
      do rowQ = 1, h%row%vec%template%noQuantities
        do rowI = 1, h%row%vec%quantities(rowQ)%template%noInstances
          rowBlock = FindBlock ( h%row, rowQ, rowI )
          do colQ1 = 1, h%col%vec%template%noQuantities
            do colI1 = 1, h%col%vec%quantities(colQ1)%template%noInstances
              colBlock1 = FindBlock ( h%col, colQ1, colI1 )
              do colQ2 = 1, h%col%vec%template%noQuantities
                do colI2 = 1, h%col%vec%quantities(colQ2)%template%noInstances
                  colBlock2 = FindBlock ( h%col, colQ2, colI2 )
                  if ( h%block ( rowBlock, colBlock1, colBlock2 ) % kind /= H_Absent ) then
                    rowPack ( rowQ ) = .true.
                    colPack ( colQ1 ) = .true.
                    colPack ( colQ2 ) = .true.
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    ! Go back and set to true any we know we want to keep
    if ( associated ( dontPack ) ) then
      do colQ = 1, j%col%vec%template%noQuantities
        if ( any ( dontPack == j%col%vec%template%quantities(colQ) ) ) &
          & colPack ( colQ ) = .true.
      end do
      do rowQ = 1, j%row%vec%template%noQuantities
        if ( any ( dontPack == j%row%vec%template%quantities(rowQ) ) ) &
          & rowPack ( rowQ ) = .true.
      end do
    end if

    ! Now work out the block mappings
    rowBlockFlag = rowPack ( j%row%quant )
    do rowBlock = 1, j%row%nb
      rowBlockMap ( rowBlock ) = count ( rowBlockFlag ( 1:rowBlock ) )
    end do
    colBlockFlag = colPack ( j%col%quant )
    do colBlock = 1, j%col%nb
      colBlockMap ( colBlock ) = count ( colBlockFlag ( 1:colBlock ) )
    end do

  end subroutine MakeMatrixPackMap

  ! --------------------------------------- PopulateL2PCBinByName ---
  subroutine PopulateL2PCBinByName ( name, IgnoreHessian )
    integer, intent(in) :: NAME         ! Name index
    logical, intent(in), optional :: IgnoreHessian
    integer :: INDEX
    ! Executable code
    index = FindFirst ( l2pcDatabase%name, name )
    if ( index /= 0 ) call PopulateL2PCBin ( index, ignoreHessian )
  end subroutine PopulateL2PCBinByName

  ! --------------------------------------- Populate L2PCBin --------
  subroutine PopulateL2PCBin ( bin, IgnoreHessian )
    use HDF5, only: hsize_t, h5gclose_f, h5gopen_f, h5gget_obj_info_idx_f, &
        & h5gn_members_f
    use MLSHDF5, only: getHDF5Attribute, getHDF5DSDims, loadFromHDF5DS

    integer, intent(in) :: BIN ! The bin index to populate
    logical, intent(in), optional :: IgnoreHessian

    ! Local variables
    type(L2pc_t), pointer :: L2PC  ! This l2pc
    type(L2PCInfo_T), pointer :: INFO ! Info for this l2pc
    type(HessianElement_T), pointer :: H0 
    type(MatrixElement_T), pointer :: M0 
    integer :: BLOCKCOL  ! Loop counter
    character ( len=64 ) :: NAME ! Name of a block
    integer :: BLOCKID   ! Group ID for a block
    integer :: BLOCKLAYER  ! Loop counter
    integer :: BLOCKROW  ! Loop counter
    integer :: BLOCKSID  ! Group ID for all the (non-)Hessian blocks
    integer, dimension(3) :: DIMS
    integer(kind=hSize_t), dimension(3) :: HDIMS ! dimensions for l2pc%h%values
    integer :: KIND      ! Kind for this block
    logical :: NoHessian ! Don't try to get the Hessian
    integer :: NOVALUES  ! Number of values for this block
    integer :: STATUS    ! Flag from HDF5
    integer :: i, nmembers, objtype
    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    verboser = switchDetail( switches, 'l2pc' ) > 0
    l2pc => l2pcDatabase ( bin )
    if ( verbose ) then
      call output( 'Populating l2pc', advance='yes' )
      call outputNamedValue( 'size(l2pcDatabase)', size(l2pcDatabase) )
      call outputNamedValue( 'size(l2pcInfo)', size(l2pcInfo) )
      call outputNamedValue( 'm_unknown', m_unknown )
      call dump( l2pc%j%block%kind, 'l2pc%j%block%kind' )
      call outputNamedValue( 'any ( l2pc%j%block%kind == m_unknown )', any ( l2pc%j%block%kind == m_unknown ) )
    endif
    if ( .not. any ( l2pc%j%block%kind == m_unknown ) ) return
    info => l2pcInfo ( bin )

    call h5gn_members_f( info%binID, 'Blocks', nmembers, status )
    if ( verbose ) then
      call outputNamedValue ( 'nmembers under ' // 'Blocks', nmembers )
      do i=0, nmembers-1
        call h5gget_obj_info_idx_f( info%binID, 'Blocks', i, &
                                           name, objtype, status )
        call outputNamedValue ( 'member name: ', trim(name) )
      enddo
    endif

    call h5gOpen_f ( info%binId, 'Blocks', blocksId, status )
    if ( status /= 0 ) then
      call outputNamedValue( 'binID', info%binID )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open Blocks group for l2pc matrix block '//&
        & trim(info%matrixName) )
    endif
    do blockRow = 1, l2pc%j%row%NB
      do blockCol = 1, l2pc%j%col%NB
        ! Skip blocks we know about or are absent
        m0 => l2pc%j%block ( blockRow, blockCol )
        if ( m0%kind /= m_unknown ) cycle
        ! Access this block
        ! write ( name, * ) 'Block', blockRow, blockCol
        call CreateBlockName ( blockRow, blockCol, 0, name )
        ! This doesn't work any more (Why not?)
        ! call h5gOpen_f ( info%blocksId, trim(name), blockId, status )
        info%blocksId = blocksID ! When did these diverge?
        call h5gOpen_f ( blocksId, trim(name), blockId, status )
        if ( status /= 0 ) then
          call outputNamedValue( 'Block name', trim(name) )
          call outputNamedValue( 'Blocks ID', info%blocksId )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to open group for l2pc matrix block '//trim(name) )
        elseif ( verboser ) then
          call outputNamedValue( 'Block name', trim(name) )
          call outputNamedValue( 'Blocks ID', info%blocksId )
          call outputNamedValue( 'Block ID', blockId )
        endif
        ! Get kind of block
        call GetHDF5Attribute ( blockID, 'kind', kind )

        ! Read the block
        if ( kind == m_banded .or. kind == m_column_sparse ) then
          call GetHDF5Attribute ( blockID, 'noValues', noValues )
          call CreateBlock ( l2pc%j, blockRow, blockCol, kind, noValues )
          m0 => l2pc%j%block ( blockRow, blockCol )
          call LoadFromHDF5DS ( blockId, 'r1', m0%r1 )
          call LoadFromHDF5DS ( blockId, 'r2', m0%r2 )
        else
          call CreateBlock ( l2pc%j, blockRow, blockCol, kind )
          m0 => l2pc%j%block ( blockRow, blockCol )
        end if
        if ( kind /= m_absent ) &
          call LoadFromHDF5DS ( blockID, 'values', m0%values )
        call h5gClose_f ( blockId, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close group for l2pc matrix block '//trim(name) )
      end do
    end do
    call h5gClose_f ( blocksId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close Blocks group for l2pc matrix '//trim(info%matrixName) )

    noHessian = .false. ! Assume we want the Hessian
    if ( present(ignoreHessian) ) noHessian = ignoreHessian
    if ( l2pc%goth .and. .not. noHessian ) then
      call h5gOpen_f ( info%binId, 'HessianBlocks', blocksId, status )
      if ( status /= 0 ) then
        call outputNamedValue( 'binID', info%binID )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open HessianBlocks group for l2pc matrix block '//&
          & trim(info%matrixName) )
      endif
      ! Now loop over the blocks and read them.
      do blockRow = 1, l2pc%h%row%NB
        do blockCol = 1, l2pc%h%col%NB
          do blockLayer = 1, l2pc%h%col%NB
            ! Skip blocks we know about or are absent
            h0 => l2pc%h%block ( blockRow, blockCol, blockLayer )
            if ( h0%kind /= h_unknown ) cycle
            call CreateBlockName ( blockRow, blockCol, blockLayer, name )
            call h5gOpen_f ( blocksId, trim(name), blockId, status )
            if ( status /= 0 ) then
              call outputNamedValue( 'Block name', trim(name) )
              call outputNamedValue( 'Blocks ID', blocksId )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to open group for l2pc Hessian block '//trim(name) )
            elseif ( verboser ) then
              call outputNamedValue( 'Block name', trim(name) )
              call outputNamedValue( 'Block ID', blockId )
              call outputNamedValue( 'HessianBlocks ID', blocksId )
            endif
            call GetHDF5Attribute ( blockID, 'kind', kind )
            if ( DEEBUG .or. verboser ) &
              & call outputNamedValue( 'kind', kind )
            select case (kind)
            case (h_absent)
              call CreateBlock ( l2pc%h, blockRow, blockCol, blockLayer, kind )
              h0 => l2pc%h%block ( blockRow, blockCol, blockLayer )
              ! call LoadFromHDF5DS ( blockId, 'i', h0%values )
            case (h_full)
              call CreateBlock ( l2pc%h, blockRow, blockCol, blockLayer, kind )
              h0 => l2pc%h%block ( blockRow, blockCol, blockLayer )
              if ( any(shape(h0%values) == 0) ) then
                call GetHDF5DSDims ( blockId, 'values', HDIMS )
                dims = hdims
                call DeAllocate_test ( h0%values, &
                  & 'h0%values', moduleName // 'PopulateOneL2PC' )
                call Allocate_test ( h0%values, dims(1), dims(2), dims(3),&
                  & 'h0%values', moduleName // 'PopulateOneL2PC' )
              endif
              call LoadFromHDF5DS ( blockId, 'values', h0%values )
            case (h_sparse)
              call GetHDF5Attribute ( blockID, 'noValues', noValues )
              if ( DEEBUG ) call outputNamedValue( 'noValues', noValues )
              call CreateBlock ( l2pc%h, blockRow, blockCol, blockLayer, kind, noValues )
              h0 => l2pc%h%block ( blockRow, blockCol, blockLayer )
              call LoadFromHDF5DS ( blockId, 'i', h0%tuples(1:h0%tuplesFilled)%i )
              call LoadFromHDF5DS ( blockId, 'j', h0%tuples(1:h0%tuplesFilled)%j )
              call LoadFromHDF5DS ( blockId, 'k', h0%tuples(1:h0%tuplesFilled)%k )
              call LoadFromHDF5DS ( blockId, 'h', h0%tuples(1:h0%tuplesFilled)%h )
            case default
              call outputNamedValue( 'h_full', h_full )
              call outputNamedValue( 'h_sparse', h_sparse )
              call outputNamedValue( 'h_unknown', h_unknown )
              call outputNamedValue( 'kind', kind )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unexpected kind for l2pc Hessian block '//trim(name) )
            end select
            call h5gClose_f ( blockId, status )
            if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to close group for l2pc Hessian matrix block '//trim(name) )
          end do
        end do
      end do
      call h5gClose_f ( blocksId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close HessianBlocks group for l2pc matrix '//trim(info%matrixName) )
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Fully populating L2PC Bin when we have a Hessian matrix is untested' )
      !!! stop !!! CODE HERE
    endif

    if ( switchdetail ( switches, 'spa' ) > -1 ) call dump_struct ( l2pc%j, 'Populated l2pc bin' )

  end subroutine PopulateL2PCBin

  ! --------------------------------------- ReadCompleteHDF5L2PCFile -------
  subroutine ReadCompleteHDF5L2PCFile ( MLSFile, Where, Shallow )
    use HDF, only: DFAcc_rdonly
    use HDF5, only: H5GN_Members_f
    use intrinsic, only: l_hdf
    use MLSFiles, only: HDFVersion_5, MLS_Openfile
    use Trace_M, only: trace_begin, trace_end
    use Toggles, only: toggle, gen
    type (MLSFile_T), intent(inout)   :: MLSFile
    integer, intent(in) :: Where ! In the L2CF tree, for tracing
    logical, optional   :: Shallow
    ! character (len=*), intent(in) :: FILENAME

    ! Local variables
    integer :: BIN             ! Loop counter
    integer :: DUMMY           ! Ignored return from AddToDatabase
    integer :: FILEID          ! From hdf5
    type (L2PCInfo_T) :: INFO  ! Info for one bin
    type (L2pc_t) :: L2PC      ! The l2pc read from one bin
    integer :: Me = -1         ! String index for trace
    logical :: MYSHALLOW       ! Value of shallow
    integer :: NOBINS          ! Number of bins
    integer :: STATUS          ! Flag from HDF5

    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    if ( verbose ) &
      & call output ( 'Reading complete l2pc file ', advance='yes' )
    myShallow = .true.
    if ( present ( shallow ) ) myShallow = shallow

    call trace_begin ( me, "ReadCompleteHDF5L2PCFile", where, cond=toggle (gen) )
    MLSFile%content = 'l2pc'
    MLSFile%access = DFACC_RDONLY
    MLSFile%HDFVersion = HDFVERSION_5
    MLSFile%type = l_hdf
    status = 0
    if ( .not. MLSFile%stillOpen ) call mls_openFile( MLSFile, status )
    fileID = MLSFile%FileID%f_id
    ! call h5fopen_f ( MLSFile%name, H5F_ACC_RDONLY_F, fileID, status )
    ! MLSFile%FileID%f_id = fileID
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open hdf5 l2pc file for input:'//trim(MLSFile%name), &
      & MLSFile=MLSFile )

    ! Get the number of bins
    call h5gn_members_f ( fileID, '/', noBins, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get number of bins from input l2pc file:'//trim(MLSFile%name), &
      & MLSFile=MLSFile )

    if ( verbose ) then
      call output ( 'Reading l2pc ' )
      call output ( trim(MLSFile%name), advance='yes' )
      call output ( 'Number of bins: ' )
      call output ( noBins, advance='yes' )
    endif

    ! Don't forget HDF5 numbers things from zero
    do bin = 0, noBins-1
      call ReadOneHDF5L2PCRecord ( L2PC, MLSFile, bin, &
        & shallow=myShallow, info=Info )
      dummy = AddL2PCToDatabase ( l2pcDatabase, L2PC )
      dummy = AddFileIDToDatabase ( fileIDDatabase, MLSFile%fileID%f_id )
      if ( switchDetail ( switches, 'spa' ) > -1 ) call Dump_struct ( l2pc%j, 'One l2pc bin' ) 

      ! Now nullify the pointers in l2pc so we don't clobber the one we've written
      nullify ( l2pc%j%block )
      nullify ( l2pc%j%row%nelts, l2pc%j%row%inst, l2pc%j%row%quant )
      nullify ( l2pc%j%col%nelts, l2pc%j%col%inst, l2pc%j%col%quant )
      nullify ( l2pc%j%row%vec%template%quantities, l2pc%j%col%vec%template%quantities )
      nullify ( l2pc%j%row%vec%quantities, l2pc%j%col%vec%quantities )

      nullify ( l2pc%h%block )
      nullify ( l2pc%h%row%nelts, l2pc%h%row%inst, l2pc%h%row%quant )
      nullify ( l2pc%h%col%nelts, l2pc%h%col%inst, l2pc%h%col%quant )
      nullify ( l2pc%h%row%vec%template%quantities, l2pc%h%col%vec%template%quantities )
      nullify ( l2pc%h%row%vec%quantities, l2pc%h%col%vec%quantities )
      dummy = AddL2PCInfoToDatabase ( l2pcInfo, Info )
    end do

    ! Don't close the file, we're keeping it open to read blocks from it later
    call trace_end ( "ReadCompleteHDF5L2PCFile", cond=toggle (gen) )
  end subroutine ReadCompleteHDF5L2PCFile

  ! --------------------------------------- ReadOneHDF5L2PCRecord ------------
  subroutine ReadOneHDF5L2PCRecord ( l2pc, MLSFile, l2pcIndex, shallow, info )
    use HDF5, only: H5GClose_f, h5gopen_f, h5gget_obj_info_idx_f, &
      & h5gn_members_f
    use MLSHDF5, only: GetHDF5Attribute, loadFromHDF5DS
    type ( L2pc_t ), intent(out), target :: L2PC
    type (MLSFile_T), intent(inout)   :: MLSFile
    ! integer, intent(in) :: FILEID       ! HDF5 ID of input file
    integer, intent(in) :: L2PCINDEX        ! Index of l2pc entry to read
    logical, optional, intent(in) :: SHALLOW ! Don't read blocks
    type(L2PCInfo_T), intent(out), optional :: INFO ! Information output

    ! Local variables
    integer :: BLOCKID                  ! Id of block group (one block at a time)
    integer :: BLOCKSID                 ! Id of blocks group
    integer :: HBLOCKSID                ! Id of HessianBlocks group
    integer :: I, J, K                  ! Block loop indices
    integer :: KIND                     ! Kind of block (absent, etc.)
    integer :: MATRIXID                 ! HDF5 for matrix group
    integer :: NMEMBERS
    integer :: NOVALUES                 ! For banded or column sparse cases
    integer :: OBJTYPE                  ! From HDF5
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGINDEX              ! Index of string
    integer :: XSTAR                    ! Linearisation state vector index
    integer :: YSTAR                    ! Radiances for xStar vector index

    logical :: COLINSTANCEFIRST         ! Flag for matrix
    logical :: ROWINSTANCEFIRST         ! Flag for matrix
    logical :: MYSHALLOW                ! Value of shallow

    type (MatrixElement_T), pointer :: M0 ! A Matrix block within kStar
    type (HessianElement_T), pointer :: H0 ! A Hessian block within kStar
    character ( len=64 ) :: MATRIXNAME  ! Name for matrix
    character ( len=64 ) :: NAME        ! Name for block group


    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    verboser = switchDetail( switches, 'l2pc' ) > 0
    myShallow = .false.
    if ( present ( shallow ) ) myShallow = shallow

    if ( verbose ) &
      & call output ( 'Reading bin from l2pc file', advance='yes' )

    call h5gGet_obj_info_idx_f ( MLSFile%fileID%f_id, '/', l2pcIndex, matrixName, &
      & objType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get information on matrix in input l2pc file', &
      & MLSFile=MLSFile )
    if ( verbose ) &
      & call outputNamedValue ( 'Reading matrix', matrixName )
    call h5gOpen_f ( MLSFile%fileID%f_id, trim(matrixName), matrixId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open matrix in input l2pc file' , &
      & MLSFile=MLSFile )
    ! Read the row and column vectors
    MLSFile%fileID%grp_id = matrixId
    call ReadOneVectorFromHDF5 ( MLSFile, 'Columns', xStar )
    call ReadOneVectorFromHDF5 ( MLSFile, 'Rows', yStar )
    if ( verbose ) then
      call output( 'x* vector', advance='yes' )
      call dump( L2PCVS(xStar) )
      call output( 'y* vector', advance='yes' )
      call dump( L2PCVS(yStar) )
    endif

    ! Get the instance first information
    call h5gOpen_f ( matrixID, 'Blocks', blocksID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to access Blocks group of input l2pc matrix' , &
      & MLSFile=MLSFile )
    ! This next is inelegant: grp_id should be a linked list
    MLSFile%fileID%sd_id = blocksID
    call GetHDF5Attribute ( MLSFile, 'rowInstanceFirst', rowInstanceFirst )
    call GetHDF5Attribute ( MLSFile, 'colInstanceFirst', colInstanceFirst )

    stringIndex = GetStringIndexFromString ( "'"//trim(matrixName)//"'" )

    ! Create the matrix
    call CreateEmptyMatrix ( l2pc%j, stringIndex, l2pcVs(yStar), l2pcVs(xStar), &
      & row_quan_first = .not. rowInstanceFirst,&
      & col_quan_first = .not. colInstanceFirst )

    ! Must *not* forget to name l2pc 
    ! (l2pc%name will itself never be used beyond this point)
    l2pc%name = l2pc%j%name

    ! Fill up the information
    if ( present ( info ) ) then
      info%fileID = MLSFile%fileID%f_id
      info%binID = matrixID
      info%blocksID = blocksID
      info%matrixName = matrixName
    end if

    if ( verboser ) then
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          call CreateBlockName ( i, j, 0, name )
          call h5gOpen_f ( blocksId, trim(name), blockId, status )
          if ( status /= 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open group for l2pc matrix block '//trim(name) , &
              & MLSFile=MLSFile )
          elseif ( verboser ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call outputNamedValue( 'Block ID', blockId )
          endif
          call h5gClose_f ( blockId, status )
        enddo
      enddo
    endif
    ! Loop over blocks and read them
    if ( .not. myShallow ) then
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          ! Access this block
          call CreateBlockName ( i, j, 0, name )
          call h5gOpen_f ( blocksId, trim(name), blockId, status )
          if ( status /= 0 ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open group for l2pc matrix block '//trim(name) , &
              & MLSFile=MLSFile )
          elseif ( verboser ) then
            call outputNamedValue( 'Block name', trim(name) )
            call outputNamedValue( 'Blocks ID', blocksId )
            call outputNamedValue( 'Block ID', blockId )
          endif
          ! Could check it's the block we're expecting but I think I'll be lazy
          ! This next is inelegant: grp_id should be a linked list
          MLSFile%fileID%sd_id = blockID
          call GetHDF5Attribute ( MLSFile, 'kind', kind )
          if ( kind == m_banded .or. kind == m_column_sparse ) then
            call GetHDF5Attribute ( MLSFile, 'noValues', noValues )
            call CreateBlock ( l2pc%j, i, j, kind, noValues )
            m0 => l2pc%j%block ( i, j )
            call LoadFromHDF5DS ( MLSFile, 'r1', m0%r1 )
            call LoadFromHDF5DS ( MLSFile, 'r2', m0%r2 )
          else
            call CreateBlock ( l2pc%j, i, j, kind )
            m0 => l2pc%j%block ( i, j )
          end if
          if ( kind /= m_absent ) &
            call LoadFromHDF5DS ( MLSFile, 'values', m0%values )
          call h5gClose_f ( blockId, status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to close group for input l2pc matrix block '//trim(name) , &
            & MLSFile=MLSFile )
        end do
      end do
    else
      ! Otherwise, flag the whole matrix as unknown
      do i = 1, l2pc%j%row%NB
        do j = 1, l2pc%j%col%NB
          call CreateBlock ( l2pc%j, i, j, m_unknown )
        end do
      end do
    end if

    call h5gClose_f ( blocksId, status )
    ! Look for any Hessian blocks
    ! We begin by assuming there are none
    call h5gn_members_f( MLSFile%fileID%f_id, trim(matrixName), nmembers, status )
    if ( verbose ) then
      call outputNamedValue ( 'nmembers under ' // trim(matrixName), nmembers )
      do i=0, nmembers-1
        call h5gget_obj_info_idx_f( MLSFile%fileID%f_id, trim(matrixName), i, &
                                           name, objtype, status )
        call outputNamedValue ( 'member name: ', trim(name) )
      enddo
    endif
    status = 1
    if ( nmembers > 3 ) call h5gOpen_f ( matrixID, 'HessianBlocks', hBlocksID, status )
    l2pc%goth = status == 0 .and. nmembers > 3 ! Got Hessian?
    if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call outputNamedValue( 'got hessian?', l2pc%goth )
    if ( l2pc%goth ) then
      ! This next is inelegant: grp_id should be a linked list
      MLSFile%fileID%sd_id = hblocksID

      l2pc%h = CreateEmptyHessian ( stringIndex, l2pcVs(yStar), l2pcVs(xStar) )
      if ( switchDetail( switches, 'hess' ) > -1 ) &
                & call output( 'Creating empty hessian', advance='yes' )
      if ( present ( info ) ) info%hBlocksID = hBlocksID

      ! Loop over blocks and read them
      if ( .not. myShallow ) then
        if ( switchDetail( switches, 'hess' ) > -1 ) then
          call output( 'trying to read this hessian block', advance='yes' )
          call outputNamedValue( 'l2pc%h%row%NB', l2pc%h%row%NB )
          call outputNamedValue( 'l2pc%h%col%NB', l2pc%h%col%NB )
          call outputNamedValue( 'hblocksId', hblocksId )
        endif
        do i = 1, l2pc%h%row%NB
          do j = 1, l2pc%h%col%NB
            do k = 1, l2pc%h%col%NB
              ! Access this block
              call CreateBlockName ( i, j, k, name )
              call h5gOpen_f ( hblocksId, trim(name), blockId, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to open group for l2pc hessian block '//trim(name) , &
                & MLSFile=MLSFile )
              ! Could check it's the block we're expecting but I think I'll be lazy
              ! This next is inelegant: grp_id should be a linked list
              MLSFile%fileID%sd_id = blockID
              call GetHDF5Attribute ( MLSFile, 'kind', kind )
              if ( switchDetail( switches, 'hess' ) > 1 ) &
                & call outputNamedValue( 'kind', kind )

              select case ( kind )
              case ( h_absent )
                call CreateBlock ( l2pc%h, i, j, k, kind, noValues )
                h0 => l2pc%h%block ( i, j, k )
              case ( h_sparse )
                call GetHDF5Attribute ( MLSFile, 'noValues', noValues )
                call CreateBlock ( l2pc%h, i, j, k, kind, noValues )
                h0 => l2pc%h%block ( i, j, k )
                call LoadFromHDF5DS ( MLSFile, 'i', h0%tuples%i )
                call LoadFromHDF5DS ( MLSFile, 'j', h0%tuples%j )
                call LoadFromHDF5DS ( MLSFile, 'k', h0%tuples%k )
                call LoadFromHDF5DS ( MLSFile, 'h', h0%tuples%h )
              case ( h_full )
                call CreateBlock ( l2pc%h, i, j, k, kind )
                h0 => l2pc%h%block ( i, j, k )                
                call LoadFromHDF5DS ( MLSFile, 'values', h0%values )
              case default
                 call MLSMessage ( MLSMSG_Warning, ModuleName, &
                & 'Read unexpected Hessian kind '//trim(name) )
              end select
              call h5gClose_f ( blockId, status )
              if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Unable to close group for input l2pc Hessian block '//trim(name) , &
                & MLSFile=MLSFile )
            end do
          end do
        end do
      else
        if ( switchDetail( switches, 'hess' ) > -1 ) &
                  & call output( 'flagging hessian as unknown', advance='yes' )
        ! Otherwise, flag the whole hessian as unknown
        do i = 1, l2pc%h%row%NB
          do j = 1, l2pc%h%col%NB
            do k = 1, l2pc%h%col%NB
              call CreateBlock ( l2pc%h, i, j, k, h_unknown )
            end do
          end do
        end do
      end if
    end if                              ! Looking for Hessians

    ! finish up, though if in shallow mode, then keep the groups open
    if ( .not. myShallow ) then
      if ( l2pc%goth ) then
        call h5gClose_f ( hblocksID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close Blocks group for input l2pc' , &
          & MLSFile=MLSFile )
      endif
      call h5gClose_f ( matrixID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close matrix group for input l2pc' , &
        & MLSFile=MLSFile )
    endif
    i = switchDetail( switches, 'l2pc' )
    if ( i >= 0 ) then
      call dump ( l2pc, details=i-1 )
      call output ( 'Done reading bin from l2pc file', advance='yes' )
    end if

  end subroutine ReadOneHDF5L2PCRecord

  ! --------------------------------------- ReadOneVectorFromHDF5 ------
  subroutine ReadOneVectorFromHDF5 ( MLSFile, name, vector )
    use HDF5, only: H5GClose_f, h5gopen_f, h5gget_obj_info_idx_f
    use MLSHDF5, only: GETHDF5Attribute, isHDF5AttributePresent, &
      & isHDF5DSPresent, loadFromHDF5DS
    use MLSSignals_m, only: radiometers, signals
    ! Read a vector from an l2pc HDF5 and adds it to internal databases.
    ! Dummy arguments
    type (MLSFile_T) :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of vector
    integer, intent(out) :: VECTOR      ! Index of vector read in L2PCVectors

    ! Local variables
    integer :: FREQUENCYCOORDINATE      ! Enumeration
    integer :: I                        ! Index into various arrays
    integer :: MOLECULE                 ! Enumeration
    integer :: NAMEINDEX                ! Quantity name
    integer :: NOCHANS                  ! Dimension
    integer :: NOINSTANCES              ! Dimension
    integer :: NOQUANTITIES             ! Number of quantities in vector
    integer :: NOSURFS                  ! Dimension
    integer :: OBJTYPE                  ! Irrelevant argument to HDF5
    integer :: QID                      ! HDF5 ID of quantity group
    integer :: QT0                      ! Offset into l2pc quantity templates database
    integer :: QTINDEXOFFSET            ! First free index in quantity template database
    integer :: QUANTITY                 ! Loop inductor
    integer :: QUANTITYTYPE             ! Enumerated
    integer :: RADIOMETER               ! Enumerated
    integer :: SIDEBAND                 ! Sideband -1,0,1
    integer :: SIGNAL                   ! Index
    integer :: STATUS                   ! Flag from HDF
    integer :: STRINGINDEX              ! A string index
    integer :: UNIT                     ! Units for quantity
    integer :: VERTICALCOORDINATE       ! Enumeratire
    integer :: VID                      ! HDF5 ID of vector group
    integer :: VTINDEX                  ! Index of vector template
    character (len=64) :: THISNAME      ! Name of this quantity
    character (len=64) :: WORD          ! General string

    logical :: COHERENT                 ! Flag
    logical :: DEBUG                    ! Flag
    logical :: STACKED                  ! Flag
    logical :: LOGBASIS                 ! Flag

    integer, dimension(:), pointer :: SIGINDS  ! Index into signals database
    integer, dimension(:), pointer :: QTINDS   ! Quantity indices

    character (len=64), pointer, dimension(:) :: QUANTITYNAMES ! Names of quantities
    type ( QuantityTemplate_T), pointer :: QT    ! Template for the quantity
    type ( VectorTemplate_T) :: VT    ! Template for the vector
    type ( Vector_T ) :: V              ! The vector

    ! Executable code
    verbose = switchDetail( switches, 'l2pc') > -1
    debug = switchDetail( switches, 'debug') > -1
    nullify ( sigInds, qtInds )

    if ( verbose ) then
      call output ( 'Reading ' )
      call output ( trim(name) )
      call output ( ' vector from l2pc', advance='yes' )
    end if
    ! Open the vector group
    call h5gOpen_f ( MLSFile%fileID%grp_id, name, vId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group for vector '//trim(name), &
      & MLSFile=MLSFile )

    ! Get the number of quantities
    MLSFile%fileID%sd_id = vID
    call GetHDF5Attribute ( MLSFile, 'noQuantities', noQuantities )
    call allocate_test ( qtInds, noQuantities, 'qtInds', ModuleName )
    ! Work out the order of the quantities
    nullify ( quantityNames )
    call Allocate_test ( quantityNames, noQuantities, 'quantityNames', ModuleName )
    do quantity = 1, noQuantities
      if ( verbose) &
        & call output ( quantity, before='Identifying quantity ', advance='yes' )
      call h5gget_obj_info_idx_f ( MLSFile%fileID%grp_id, name, quantity-1, thisName, &
        & objType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to access a quantity within '//trim(name) , &
      & MLSFile=MLSFile )
      call h5gOpen_f ( vId, trim(thisName), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
      MLSFile%fileID%sd_id = qID
      call GetHDF5Attribute ( MLSFile, 'index', i )
      quantityNames ( i ) = thisName
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    qtIndexOffset = InflateQuantityTemplateDatabase ( l2pcQTs, noQuantities )
    if ( debug ) then
      call outputNamedValue( 'noQuantities', noQuantities )
      call outputNamedValue( 'qt0', qt0 )
    endif
    qt0 = size ( l2pcQTs ) - noQuantities

    ! Now go through quantities in order
    do quantity = 1, noQuantities
      qt => l2pcQTs ( qtIndexOffset + quantity - 1 )
      if ( verbose ) then
        call output ( quantity, before='Reading quantity ' )
        call output ( ': ' )
        call output ( trim(quantityNames(quantity)), advance='yes' )
      end if

      ! Point to this quantity
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name) , &
        & MLSFile=MLSFile )

      ! Store the name
      stringIndex = GetStringIndexFromString ( quantityNames(quantity) )
      nameIndex = stringIndex

      ! Get the quantity type
      MLSFile%fileID%sd_id = qID
      call GetHDF5Attribute ( MLSFile, 'type', word )
      quantityType = GetLitIndexFromString ( word )

      ! Get other info as appropriate
      signal = 0
      sideband = 0
      molecule = 0
      radiometer = 0
      frequencyCoordinate = l_none
      unit = phyq_dimensionless
      select case ( quantityType )
      case ( l_rhi )
        word = 'H2O'
        molecule = GetLitIndexFromString ( word )
        unit = phyq_pctrhi
      case ( l_columnAbundance )
        call GetHDF5Attribute ( MLSFile, 'molecule', word )
        molecule = GetLitIndexFromString ( word )
        unit = phyq_colmabundance
      case ( l_isotopeRatio )
        call GetHDF5Attribute ( MLSFile, 'molecule', word )
        molecule = GetLitIndexFromString ( word )
      case ( l_vmr )
        call GetHDF5Attribute ( MLSFile, 'molecule', word )
        molecule = GetLitIndexFromString ( word )
        if ( isExtinction(molecule) ) then
          if ( IsHDF5AttributePresent ( MLSFile, 'radiometer' ) ) then
            call GetHDF5Attribute ( MLSFile, 'radiometer', word )
            stringIndex = GetStringIndexFromString ( word )
            radiometer = FindFirst ( Radiometers%prefix, stringIndex )
          else
            radiometer = 0
          endif
          frequencyCoordinate = l_intermediateFrequency
        else
          frequencyCoordinate = l_none
        end if
        unit = phyq_vmr
      case ( l_radiance, l_TScat, l_calsidebandfraction, l_limbsidebandfraction, &
        & l_opticalDepth, l_elevOffset)
        call GetHDF5Attribute ( MLSFile, 'signal', word )
        if ( quantityType /= l_TScat ) then
          call Parse_Signal ( word, sigInds, sideband=sideband )
        else ! TScat
          call Parse_Signal ( word, sigInds, sideband=sideband, channels=qt%channels )
          call allocate_test ( qt%chanInds, count(qt%channels), 'qt%chanInds', moduleName )
          qt%chanInds = pack( (/ ( i, i=lbound(qt%channels,1), ubound(qt%channels,1) ) /), &
                            & qt%channels ) ! Indices of true channels
          call allocate_test ( qt%frequencies, size(qt%chanInds), 'qt%frequencies', moduleName )
          qt%frequencies = signals(sigInds(1))%frequencies(qt%chanInds)
        end if
        signal = sigInds(1)
        verticalCoordinate = l_geodAltitude
        frequencyCoordinate = l_channel
        call deallocate_test(sigInds,'sigInds',ModuleName)
        unit = phyq_temperature
      case ( l_temperature )
        unit = phyq_temperature
      case default
      end select

      ! Now read the dimensions for the quantity
      call GetHDF5Attribute ( MLSFile, 'noChans', noChans )
      call GetHDF5Attribute ( MLSFile, 'noInstances', noInstances )
      call GetHDF5Attribute ( MLSFile, 'noSurfs', noSurfs )
      call GetHDF5Attribute ( MLSFile, 'verticalCoordinate', word )
      verticalCoordinate = GetLitIndexFromString ( word )
      ! Look for frequency coordinate (optional for backwards compatability with
      ! older l2pc files.
      if ( IsHDF5AttributePresent ( qID, 'frequencyCoordinate' ) ) then
        call GetHDF5Attribute ( MLSFile, 'frequencyCoordinate', word )
        frequencyCoordinate = GetLitIndexFromString ( trim(word) )
      end if
      call GetHDF5Attribute ( MLSFile, 'logBasis', logBasis )
      call GetHDF5Attribute ( MLSFile, 'coherent', coherent )
      call GetHDF5Attribute ( MLSFile, 'stacked', stacked )

      qt%name = nameIndex
      call SetupNewQuantityTemplate ( qt, &
        & noInstances=noInstances, noSurfs=noSurfs, noChans=noChans, &
        & coherent=coherent, stacked=stacked )

      qt%quantityType = quantityType
      qt%logBasis = logBasis
      qt%molecule = molecule
      qt%radiometer = radiometer
      qt%signal = signal
      qt%sideband = sideband
      qt%verticalCoordinate = verticalCoordinate
      qt%frequencyCoordinate = frequencyCoordinate
      qt%unit = unit
      if ( debug ) then
        call outputNamedValue( 'quantity', quantity )
        call outputNamedValue( 'quantityType', quantityType )
        call outputNamedValue( 'logBasis', logBasis )
        call outputNamedValue( 'molecule', molecule )
        call outputNamedValue( 'radiometer', radiometer )
        call outputNamedValue( 'qtIndexOffset', qtIndexOffset )
        call outputNamedValue( 'index', qtIndexOffset + quantity - 1 )
      endif

      ! Get the surfaces and phis
      call LoadFromHDF5DS ( MLSFile, 'surfs', qt%surfs )
      call LoadFromHDF5DS ( MLSFile, 'phi', qt%phi )
      if ( IsHDF5DSPresent ( qId, 'solarZenith' ) ) then
        call LoadFromHDF5DS ( MLSFile, 'solarZenith', qt%solarZenith )
      else
        qt%solarZenith = 0.0
      endif
      ! Try to get the frequencies
      if ( IsHDF5DSPresent ( qId, 'frequencies' ) ) then
        call Allocate_test( qt%frequencies, qt%noChans, &
          & 'qt%frequencies', ModuleName )
        call LoadFromHDF5DS ( MLSFile, 'frequencies', qt%frequencies )
      end if

      ! Now record the index for this quantity template
      qtInds ( quantity ) = quantity + qt0

      ! For the moment, close the quantity. We'll come back to it later to fill
      ! up the values for the vector
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    ! Now create a vector template with these quantities
    call ConstructVectorTemplate ( 0, l2pcQTs, qtInds, vt, forWhom=moduleName )
    vtIndex = AddVectorTemplateToDatabase ( l2pcVTs, vt )

    call deallocate_test ( qtInds, 'qtInds', ModuleName )

    ! Now create a vector for this vector template
    if ( verbose ) &
      & call output ( 'Creating vector', advance='yes' )
    v = CreateVector ( 0, l2pcVTs(vtIndex), l2pcQTs, vectorNameText='_v_'//name )
    if ( verbose ) &
      & call output ( 'Adding vector to database', advance='yes' )
    vector = AddVectorToDatabase ( l2pcVs, v )

    ! Now go through the quantities again and read the values
    do quantity = 1, noQuantities
      if ( verbose ) then
        call output ( 'Reading values for ' )
        call output ( trim(quantityNames(quantity)), advance='yes' )
      end if
      call h5gOpen_f ( vId, trim(quantityNames(quantity)), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(quantityNames(quantity))// &
        & ' in vector '//trim(name), &
        & MLSFile=MLSFile )
      MLSFile%fileID%sd_id = qID
      call LoadFromHDF5DS ( MLSFile, 'values', &
        & l2pcVs(vector)%quantities(quantity)%values )
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(thisName)//' in vector '//trim(name), &
        & MLSFile=MLSFile )
    end do

    call Deallocate_test ( quantityNames, 'quantityNames', ModuleName )

  end subroutine ReadOneVectorFromHDF5

  ! --------------------------------------- WriteVectorAsHDF5L2PC ----------
  subroutine WriteVectorAsHDF5L2PC ( location, vector, name, packInfo )
    use HDF5, only: H5GClose_f, h5gcreate_f
    use MLSHDF5, only: MakeHDF5Attribute, saveasHDF5DS
    use mlssignals_m, only: radiometers
    integer, intent(in) :: LOCATION     ! The HDF5 location for the vector
    type (Vector_T), intent(in), target :: VECTOR ! The vector to write
    character(len=*), intent(in) :: NAME ! Name for vector
    logical, dimension(:), intent(in) :: PACKINFO ! Pack vector if set

    ! Local variables
    character ( len=132 ) :: LINE       ! A line of text
    character ( len=32 ) :: QNAME       ! Name of quantity
    integer :: QID                      ! HDF5 ID for quantity group
    integer :: QINDEX                   ! Quantity index
    integer :: QUANTITY                 ! Index
    integer :: STATUS                   ! Flag from HDF5
    integer :: VID                      ! HDF5 ID for vector group
    type (QuantityTemplate_T), pointer :: QT ! The template

    ! Executable code
    ! Create a group for the vector
    call h5gCreate_f ( location, trim(name), vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for vector ' // trim(name) )
    call MakeHDF5Attribute ( vID, 'noQuantities', count(packInfo) )
    qIndex = 0
    do quantity = 1, size ( vector%quantities )
      if ( packInfo(quantity) ) then
        qIndex = qIndex + 1
        qt => vector%quantities(quantity)%template
        ! Create a group for the quantity
        call get_string ( qt%name, qName )
        call h5gCreate_f ( vID, trim(qName), qID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create group for quantity ' // trim(qName) )

        ! Write quantity type and index
        call MakeHDF5Attribute ( qID, 'index', qIndex )
        call get_string ( lit_indices(qt%quantityType), line )
        call MakeHDF5Attribute ( qID, 'type', trim(line) )
        ! Write other info as appropriate
        select case ( qt%quantityType )
        case (l_vmr, l_columnAbundance, l_rhi, l_isotopeRatio)
          call get_string ( lit_indices(qt%molecule), line )
          call MakeHDF5Attribute ( qID, 'molecule', trim(line) )
          if ( isExtinction(qt%molecule) ) then
            line = 'none'
            if ( qt%radiometer /= 0 ) &
              & call get_string ( radiometers(qt%radiometer)%prefix, line )
            call MakeHDF5Attribute ( qID, 'radiometer', trim(line) )
          end if
        case (l_radiance, l_calsidebandfraction, l_limbsidebandfraction, &
          & l_opticalDepth, l_elevOffset)
          call GetSignalName ( qt%signal, line, sideband=qt%sideband )
          call MakeHDF5Attribute ( qID, 'signal', trim(line) )
        case (l_TScat)
          call GetSignalName ( qt%signal, line, sideband=qt%sideband, &
            & otherChannels=qt%channels )
          call MakeHDF5Attribute ( qID, 'signal', trim(line) )
        end select
        ! Write out the dimensions for the quantity
        call MakeHDF5Attribute ( qID, 'noChans', qt%noChans )
        call MakeHDF5Attribute ( qID, 'noInstances', qt%noInstances )
        call MakeHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
        call get_string ( lit_indices(qt%verticalCoordinate), line )
        call MakeHDF5Attribute ( qID, 'verticalCoordinate', trim(line) )
        call get_string ( lit_indices(qt%frequencyCoordinate), line )
        call MakeHDF5Attribute ( qID, 'frequencyCoordinate', trim(line) )
        call get_string ( lit_indices(qt%unit), line )
        call MakeHDF5Attribute ( qID, 'unit', trim(line) )
        call MakeHDF5Attribute ( qID, 'logBasis', qt%logBasis )
        call MakeHDF5Attribute ( qID, 'coherent', qt%coherent )
        call MakeHDF5Attribute ( qID, 'stacked', qt%stacked )
        ! Write out important coordinates
        call SaveAsHDF5DS ( qID, 'phi', qt%phi )
        call SaveAsHDF5DS ( qID, 'surfs', qt%surfs )
        call SaveAsHDF5DS ( qID, 'solarZenith', qt%solarZenith )
        if ( associated ( qt%frequencies ) ) &
          & call SaveAsHDF5DS ( qID, 'frequencies', qt%frequencies )
        ! Write out the values
        if ( .not. associated(vector%quantities(quantity)%values) ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'values field not associated for this block'  // trim(qName) )
        else
          call SaveAsHDF5DS ( qID, 'values', vector%quantities(quantity)%values )
        endif
        ! Close the group
        call h5gClose_f ( qID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close quantity group for ' // trim(qName) )
      end if
    end do

    ! Finish off
    call h5gClose_f ( vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close vector group for ' // trim(name) )

  end subroutine WriteVectorAsHDF5L2PC

  subroutine CreateBlockName ( row, column1, column2, name )
    ! Because the Intel compiler produces a different result from
    ! write(str, *) int
    ! Args
    integer, intent(in)           :: row, column1, column2
    character(len=*), intent(out) :: name
    ! Internal variables
    character(len=16), dimension(3) :: strs
    ! Executable
    call writeIntsToChars( (/row, column1, column2/), strs )
    name = ' Block' // ' ' // trim(adjustl(strs(1))) // ' ' // trim(adjustl(strs(2)))
    if ( column2 > 0 ) name = trim ( name ) // ' ' // trim(adjustl(strs(3)))
  end subroutine CreateBlockName

  ! ----------------- getNBinsInFile
  ! find indices of the entries in the fileIDDatabase for a given L2PCFile
  subroutine getNBinsInFile( L2PCFile, indices, nBins )
    ! Args:
    type (MLSFile_T), intent(in)       :: L2PCFile
    integer, dimension(:), intent(out) :: indices
    integer, intent(out)               :: nBins
    ! Internal variables
    integer :: i
    ! Executable
    indices = 0
    nBins = 0
    if ( .not. associated(L2PCDataBase) ) return
    call outputNamedValue( 'size(L2PCDataBase)', size(L2PCDataBase) )
    do i=1, size(L2PCDataBase)
      if ( L2PCFile%FileID%f_id /= fileIDDataBase(i) ) cycle
      nBins = nBins + 1
      if ( nBins <= size(indices) ) indices(nBins) = i
    enddo
  end subroutine getNBinsInFile

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2PC_m

! $Log$
! Revision 2.133  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.132  2015/07/14 23:23:52  pwagner
! Use more robust isStringInTable
!
! Revision 2.131  2015/06/30 18:38:12  pwagner
! Should not crash just because some string not in table
!
! Revision 2.130  2015/04/28 23:59:52  pwagner
! Diffs made more manageable
!
! Revision 2.129  2015/03/28 01:10:03  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.128  2014/09/04 23:48:52  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.127  2014/01/09 18:15:42  pwagner
! Fixed bug in preventing crashes when qt%radiometer is 0
!
! Revision 2.126  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.125  2014/01/08 18:31:48  pwagner
! Fixed crash on gettting extinction radiometer prefix
!
! Revision 2.124  2013/08/30 03:56:01  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.123  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.122  2013/06/12 02:10:47  vsnyder
! Cruft removal
!
! Revision 2.121  2012/02/17 00:21:05  pwagner
! May dump layout of l2pc only with 'L' option
!
! Revision 2.120  2012/02/13 23:27:27  pwagner
! Fixed more bugs in read/write of quantities in xstar,ystar
!
! Revision 2.119  2012/02/02 01:17:03  pwagner
! Added LoadHessian; fixed errors that cause crashes
!
! Revision 2.118  2011/11/11 00:32:29  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.117  2011/11/09 17:41:02  pwagner
! Fixed bug added with last commit
!
! Revision 2.116  2011/11/04 23:38:47  pwagner
! Made consistent with MLSFile interface to MLSHDF5
!
! Revision 2.115  2011/08/20 02:04:03  vsnyder
! Remove unused use name
!
! Revision 2.114  2011/08/09 23:17:43  pwagner
! Tries to flush hessian bins, too
!
! Revision 2.113  2011/07/28 18:03:37  pwagner
! Fixed bug in FlushL2PCBins added in rev. 2.88
!
! Revision 2.112  2011/05/05 15:20:04  pwagner
! Fixed bugs in writing, reading radiometer attribute for extinctionv2
!
! Revision 2.111  2011/04/15 00:20:32  pwagner
! Protect against unlikely case where there is no Jacobian to output to file
!
! Revision 2.110  2011/02/18 17:56:16  pwagner
! Prevented crashes when run w/o l2cf
!
! Revision 2.109  2010/11/25 01:19:27  pwagner
! Fixed bugs in PopulateL2PCBin with Hessians
!
! Revision 2.108  2010/11/19 23:58:19  pwagner
! Added Diff
!
! Revision 2.107  2010/11/05 00:31:10  pwagner
! Too few args to DumpOneL2PC
!
! Revision 2.106  2010/11/03 18:30:41  pwagner
! Dumps now take option to say whether to dump hessians, matrices, and which blocks
!
! Revision 2.105  2010/09/25 01:14:34  vsnyder
! Add some TScat support, and some cannonball polishing
!
! Revision 2.104  2010/09/17 00:02:57  pwagner
! Can constrain which blocks dumped by name
!
! Revision 2.103  2010/08/31 02:05:38  vsnyder
! Don't nullufy qt%channels before allocating, to avoid leaking memory
!
! Revision 2.102  2010/08/27 21:55:35  pwagner
! WIll not try to write empty Hessian blocks
!
! Revision 2.101  2010/08/20 23:19:50  pwagner
! May specify which blocks to dump by name; dont dump empty blocks
!
! Revision 2.100  2010/08/13 22:22:59  pwagner
! details now can choose between matrices, hessians; fixed a bug in reading HessianBlocks
!
! Revision 2.99  2010/08/06 22:58:40  pwagner
! Added new 'hess' switch; using only switchDetail now
!
! Revision 2.98  2010/06/28 17:00:49  pwagner
! Fixed a few bugs
!
! Revision 2.97  2010/05/20 00:51:50  pwagner
! Can populate Hessians; untested
!
! Revision 2.96  2010/05/19 17:51:41  pwagner
! Details now used to dump or not L2PCInfo%blockID
!
! Revision 2.95  2010/05/19 00:32:06  vsnyder
! Pass IgnoreHessian into PopulateL2PCBin
!
! Revision 2.94  2010/05/13 23:45:14  pwagner
! Temporary expedients for l2pc files with Hessians; needs more code
!
! Revision 2.93  2010/04/30 22:55:32  vsnyder
! Give more meaningful names to vectors
!
! Revision 2.92  2010/04/17 01:43:35  vsnyder
! Simplify dump.  Set l2pc%goth if there is a Hessian block.
!
! Revision 2.91  2010/03/24 20:48:55  vsnyder
! Replace 'continue' by 'cycle'.  Call OptimizeBlock from WriteOneHDF5L2PC
! when writing a Hessian, to avoid trying to write empty arrays.  Add some
! comments.  Fix a wrong-subscript bug.
!
! Revision 2.90  2010/03/19 20:18:47  pwagner
! Moved DumpL2PCInfo out of public scope to appease ifort v10 which crashed the goldbrick
!
! Revision 2.89  2010/03/17 20:59:10  pwagner
! Code around bombs in DestroyL2PCInfoDatabase
!
! Revision 2.88  2010/02/25 18:04:45  pwagner
! l2pc type can now hold matrix and Hessian types
!
! Revision 2.87  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.86  2009/11/17 23:43:51  vsnyder
! Add ability to output TScat
!
! Revision 2.85  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.84  2008/12/18 21:31:18  pwagner
! May now dump an l2pc or allL2PCs (use with caution)
!
! Revision 2.83  2008/09/30 22:28:03  vsnyder
! Remove AuxGrids -- didn't need them after all
!
! Revision 2.82  2008/06/06 01:55:29  vsnyder
! Process AuxGrids field of vector quantities
!
! Revision 2.81  2008/05/03 01:49:55  vsnyder
! Repair a comment
!
! Revision 2.80  2007/10/03 23:58:46  vsnyder
! Add 'where' for tracing
!
! Revision 2.79  2007/04/26 20:30:32  pwagner
! Bugfix for way ifc writes ints to strings
!
! Revision 2.78  2006/08/05 02:11:58  vsnyder
! Add ForWhom argument to ConstructVectorTemplate
!
! Revision 2.77  2006/08/04 01:53:43  vsnyder
! Need to set the name, not quantityType before SetupNewQuantityTemplate
!
! Revision 2.76  2006/08/03 01:10:06  vsnyder
! Put l2cf names in leak track database
!
! Revision 2.75  2005/06/29 00:42:35  pwagner
! Passes MLSFiles to GetHDF5Attribute, LoadFromHDF5DS
!
! Revision 2.74  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.73  2004/06/10 00:57:47  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.72  2004/02/11 01:08:44  livesey
! Removed print statement.
!
! Revision 2.71  2004/01/30 23:27:32  livesey
! Capitalized the names of the bins on output
!
! Revision 2.70  2004/01/24 01:44:36  livesey
! Removed print statement
!
! Revision 2.69  2004/01/24 01:01:22  livesey
! Improvements to the adoption stuff
!
! Revision 2.68  2004/01/23 19:07:49  livesey
! More on adoption / loading
!
! Revision 2.67  2004/01/23 05:37:34  livesey
! Added the adoption stuff
!
! Revision 2.66  2003/09/15 17:45:03  livesey
! Added target declaration for fussy intel compiler
!
! Revision 2.65  2003/08/14 20:24:23  livesey
! Added the exact bin selector stuff
!
! Revision 2.64  2003/08/13 00:47:34  livesey
! Added the default bin selectors stuff
!
! Revision 2.63  2003/08/08 23:04:45  livesey
! Added the dontPack stuff for output.
!
! Revision 2.62  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.61  2003/05/13 04:46:24  livesey
! Renamed a routine to avoid conflict with VectorHDF5
!
! Revision 2.60  2003/02/06 01:37:12  livesey
! Set solarZenith to zero if reading from hdf5 and not present
!
! Revision 2.59  2003/02/06 01:34:21  livesey
! Added writing of solarZenith, and reading if it is present.
!
! Revision 2.58  2003/02/06 00:45:37  livesey
! Added name fragment to binSelector_t
!
! Revision 2.57  2003/02/05 21:55:45  livesey
! Bin selectors don't contain signal information anymore.
!
! Revision 2.56  2003/01/28 02:41:10  livesey
! Bug fix in writing matrix, now gets right element from database in all
! cases.
!
! Revision 2.55  2003/01/13 19:28:20  pwagner
! Moved several uses to speed Lahey buggs compiler
!
! Revision 2.54  2002/11/25 11:47:34  mjf
! Put SaveAsHDF5DS back to writing r4 reals.
!
! Revision 2.53  2002/11/22 12:47:47  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.52  2002/11/20 21:06:25  livesey
! Various bug fixes to make the packed storage work.
!
! Revision 2.51  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.50  2002/10/05 00:42:31  livesey
! Modified to use new pack/unpack literals/strings
!
! Revision 2.49  2002/10/02 23:20:19  livesey
! Added radiometer and saving as single precision
!
! Revision 2.48  2002/09/11 17:43:38  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.47  2002/08/28 21:43:38  livesey
! Cosmetic changes
!
! Revision 2.46  2002/08/28 20:42:39  livesey
! Now uses inflateQuantityTemplateDatabase, much much faster on reading
! l2pcs.
!
! Revision 2.45  2002/08/28 01:12:20  livesey
! Bug fix in frequency coordinate.
!
! Revision 2.44  2002/08/23 01:23:44  livesey
! Now optionally reads frequency coordinate and frequencies
!
! Revision 2.43  2002/08/21 23:10:11  livesey
! No longer scans each bin for blocks
!
! Revision 2.42  2002/08/20 21:01:40  livesey
! Bug fix in DestroyL2PCInfoDatabase
!
! Revision 2.41  2002/08/07 00:05:14  livesey
! Moved H5Open/Close_F into tree walker.  Made some error messages more
! informative.
!
! Revision 2.40  2002/07/25 15:53:40  mjf
! Initialised some elements of qt in ReadOneVectorFromASCII (already
! done in ReadOneVectorFromHDF5)
!
! Revision 2.39  2002/07/22 03:29:23  livesey
! Bug fix
!
! Revision 2.38  2002/07/22 03:25:23  livesey
! Bug fix and rework of loading vectors.
!
! Revision 2.37  2002/07/17 06:00:21  livesey
! Got hdf5 l2pc reading stuff working
!
! Revision 2.36  2002/07/11 22:21:00  pwagner
! These hdf5-savvy versions transferred from he5lib
!
! Revision 1.3  2002/07/09 17:38:17  livesey
! Added logBasis attribute
!
! Revision 1.2  2002/06/22 23:12:13  livesey
! Added sparsity dump
!
! Revision 1.1  2002/06/12 17:59:53  livesey
! Work on HDF5 stuff
!
! Revision 2.32  2002/05/21 01:13:24  livesey
! New file format includes name for reading into IDL
!
! Revision 2.31  2002/03/15 21:22:42  livesey
! Slight modification to binSelectors stuff
!
! Revision 2.30  2002/03/13 22:00:53  livesey
! Added output of vertical coordinate for xStar/yStar.
! Note reading routine deduces value I think.
!
! Revision 2.29  2002/03/13 01:28:48  livesey
! Commented out an error message I think I don't need
!
! Revision 2.28  2002/02/08 22:51:28  livesey
! Minor changes and tidy up
!
! Revision 2.27  2002/02/05 04:14:26  livesey
! Bug fix, still problems with packed write
!
! Revision 2.26  2002/01/23 00:50:50  pwagner
! Initialize binselectors to null
!
! Revision 2.25  2002/01/22 18:14:29  livesey
! Fixed typo.
!
! Revision 2.24  2002/01/21 23:11:49  livesey
! Completed BinSelectors support
!
! Revision 2.23  2002/01/21 21:13:42  livesey
! Added BinSelector definitions and support
!
! Revision 2.22  2002/01/18 00:34:23  livesey
! Added packed option to writeonel2pc, with supporting code.
!
! Revision 2.21  2001/06/21 22:45:43  livesey
! Found another one of those `clobbering the latest database item' gotchas in the code.
!
! Revision 2.20  2001/06/06 17:27:29  pwagner
! DEBUG parameter; checks on some more reads
!
! Revision 2.19  2001/05/23 22:00:59  livesey
! Changed counter start to avoid conflict with L2Parallel
!
! Revision 2.18  2001/05/02 20:24:03  livesey
! Removed some unused variables.
!
! Revision 2.17  2001/05/02 20:23:10  vsnyder
! Specify a name for a created vector
!
! Revision 2.16  2001/04/28 01:38:36  livesey
! Now maintains proper IDs for its own quantity and vector templates
!
! Revision 2.15  2001/04/28 01:26:50  livesey
! This one seems to work!
!
! Revision 2.14  2001/04/27 17:36:33  livesey
! An interim version
!
! Revision 2.13  2001/04/27 07:24:50  livesey
! Interim version.  Still loosing parts of kStar on add to database.
! Might this be a compiler bug?
!
! Revision 2.12  2001/04/27 07:05:28  livesey
! Not sure what happened, needed to restore the use statements.
!
! Revision 2.11  2001/04/26 23:55:17  livesey
! Interim version
!
! Revision 2.10  2001/04/26 22:32:04  vsnyder
! Alphabetize USEs and declarations
!
! Revision 2.9  2001/04/26 22:12:21  livesey
! Fixed, gets l_zeta, l_none
!
! Revision 2.8  2001/04/26 22:08:39  livesey
! Add check on vertical coordinates
!
! Revision 2.7  2001/04/26 20:02:26  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.6  2001/04/26 19:33:03  livesey
! Working version, reads and writes, (but no arithmetic :-) )
!
! Revision 2.5  2001/04/26 02:48:08  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.4  2001/04/26 00:06:33  livesey
! Many changes. Working towards a working read routine
!
! Revision 2.3  2001/04/25 20:32:42  livesey
! Interim version, tidied up write
!
! Revision 2.2  2001/04/24 20:20:48  livesey
! Word bin dropped from various places e.g. type
!
! Revision 2.1  2001/04/24 20:07:44  livesey
! Moved in from l2
!

