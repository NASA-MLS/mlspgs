! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Read_Spectroscopy_TeX_m

  implicit NONE
  private
  public :: Read_Spectroscopy_TeX

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Read_Spectroscopy_TeX ( Cross_Reference_File, Mol_Data_File, &
    & Line_Data_File, Catalog, Lines, Stat, Bands, Warn )

    ! Read the LaTeX-format spectroscopy files.
    ! String indices (catalog%species_name and catalog%line%line_name) are
    ! not filled.
    ! Signals, Sidebands, and Polarized fields of catalog%line are
    ! not filled.
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use, intrinsic :: ISO_Fortran_env, only: IOSTAT_END
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MOLECULES, only: FIRST_MOLECULE, LAST_MOLECULE, GetMoleculeIndex
    use Spectroscopy_Types, only: Catalog_t, Line_t
    character(*), intent(in) :: Cross_Reference_File ! For molecule names
    character(*), intent(in) :: Mol_Data_File        ! in LaTeX format
    character(*), intent(in) :: Line_Data_File       ! in LaTeX format
    type(catalog_t), intent(out) :: Catalog(first_molecule:last_molecule)
    type(line_t), pointer :: Lines(:)
    integer, intent(out) :: Stat ! IOSTAT from open, nonzero if one failed
    character(*), pointer, optional :: Bands(:) ! Bands for each line
    logical, intent(in), optional :: Warn ! Warn of missing molecules if present
                                          ! and true

    character, parameter :: Back = achar(92) ! Backslash

    integer :: A         ! count &'s to find bands in Lines record
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I
    integer :: IOSTAT
    integer :: L
    integer :: Mol
    logical :: MyWarn
    integer :: NLines    ! Total number of lines for all molecules
    integer :: NLinesMol ! Lines for one molecule
    integer :: N_Names
    character(len=255) :: Text
    character(24), allocatable :: TeX_name(:), XRef_name(:)
    integer :: U

    myWarn = .false.
    if ( present(warn) ) myWarn = warn

    ! Open the cross-reference file
    open ( newunit=u, file=cross_reference_file, form='formatted', &
      & status='old', iostat=iostat, iomsg=text )
    if ( iostat /= 0 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        'Opening ' // trim(cross_reference_file) // ' failed: ' // trim(text) )
      return
    end if

    ! Count the names in the cross-reference file
    n_names = 0
    ! Skip heading lines
    do
      read ( u, '(a)', end=2 ) text
      if ( text(1:9) == '----*----' ) exit
    end do
    ! Now read and count the names
    do
      read ( u, '(a)', end=2 ) text
      n_names = n_names + 1
    end do
  2 continue

    if ( n_names == 0 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, 'No cross reference names' )
      return
    end if

    ! Now read the names in the cross-reference file
    rewind ( u )
    allocate ( TeX_name(n_names), xref_name(n_names), stat=stat )
    call test_allocate ( stat, moduleName, 'TeX_name or XRef_name' )
    ! Skip heading lines
    do
      read ( u, '(a)' ) text ! A convenient unused variable
      if ( text(1:9) == '----*----' ) exit
    end do
    ! Now read the names
    do i = 1, n_names
      read ( u, '(a24,1x,a24)' ) TeX_name(i), xref_name(i)
      TeX_name(i) = adjustl(TeX_name(i))
      xref_name(i) = adjustl(xref_name(i))
    end do
    close ( u )

    ! Open the molecule data file
    open ( u, file=mol_data_file, form='formatted', status='old', &
      & iostat=iostat, iomsg=text )
    if ( iostat /= 0 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        'Opening ' // trim(mol_data_file) // ' failed: ' // trim(text) )
      return
    end if

    ! Read the molecule data file
    do
      read ( u, '(a)', end=4 ) text
      text = adjustl(text)
      if ( text(1:1) == '%' ) cycle
      i = index(text,'&')
      if ( i == 0 ) then
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Incorrect format in ' // trim(mol_data_file) // ': ' // trim(text) )
        close ( u )
        go to 9
      end if
      mol = findMolecule(text(:i-1))
      if ( mol == 0 ) then
        if ( myWarn ) call MLSMessage ( MLSMSG_Warning, moduleName, &
          & trim(text(:i-1)) // ' Not a molecule in ' // trim(mol_data_file) // ': ' // trim(text) )
        cycle
      end if
      catalog(mol)%molecule = mol
      text = text(i+1:)
      l = index(text,back//back) - 1
      if ( l < 0 ) l = len_trim(text)
      do i = 1, l
        if ( text(i:i) == '&' ) text(i:i) = ','
      end do
      ! Replace LaTeX exponents with Fortran exponents
      do
        i = index(text,'$'//back//'times 10^{')
        if ( i == 0 ) exit
        text(i:i) = 'e'
        text(i+1:) = text(i+len('.times 10^{')+1:) ! Trim out \times 10^{
        i = index(text,'}$')
        if ( i /= 0 ) text(i:) = text(i+2:)
      end do
      read ( text(:l), * ) catalog(mol)%defaultIsotopeRatio, catalog(mol)%mass, &
        & catalog(mol)%Qlog, catalog(mol)%continuum
      where ( catalog(mol)%Qlog > 0 ) &
        & catalog(mol)%Qlog = log10(catalog(mol)%Qlog)
    end do
  4 continue
    close ( u )

    ! Open the lines database file
    open ( u, file=line_data_file, form='formatted', status='old', &
      & iostat=iostat, iomsg=text )
    if ( iostat /= 0 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        'Opening ' // trim(line_data_file) // ' failed: ' // trim(text) )
      return
    end if

    ! Count the lines
    nLines = 0
    text = ''
    do
      ! Find the next molecule heading
      do while ( index(text,'multicolumn') == 0 )
        read ( u, '(a)', end=6 ) text
      end do
      call getMolName
      if ( mol < first_molecule .or. mol > last_molecule ) then
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Not a molecule in ' // trim(mol_data_file) // ': ' // trim(text) )
        close ( u )
        go to 9
      end if
      ! Read and count the lines for the molecule
      NLinesMol = 0
      do
        read ( u, '(a)', iostat=iostat ) text
        if ( iostat == iostat_end ) exit
        if ( index(text,'multicolumn') /= 0 ) exit
        nLines = nLines + 1
        nLinesMol = nLinesMol + 1
      end do
      allocate ( catalog(mol)%lines(nLinesMol), stat=stat )
      addr = 0
      if ( stat == 0 .and. nLinesMol > 0 ) &
        & addr = transfer(c_loc(catalog(mol)%lines(1)), addr)
      call test_allocate ( stat, moduleName, 'catalog(mol)%lines', (/ 1 /), &
        &  (/ nLinesMol /), storage_size(catalog) / 8, address=addr )
      if ( iostat == iostat_end ) exit
    end do
  6 rewind ( u )
    allocate ( lines(1:nLines), stat=stat )
    addr = 0
    if ( stat == 0 .and. nLines > 0 ) addr = transfer(c_loc(lines(1)), addr)
    call test_allocate ( stat, moduleName, 'Lines', (/ 1 /), (/ nLines /), &
      & storage_size(lines) / 8, address=addr )
    if ( present(bands) ) then
      allocate ( bands(1:nLines), stat=stat )
      addr = 0
      if ( stat == 0 .and. nLines > 0 ) addr = transfer(c_loc(bands(1)(1:1)), addr)
      call test_allocate ( stat, moduleName, 'Bands', (/ 1 /), (/ nLines /), &
        & storage_size(bands) / 8, address=addr )
    end if

    ! Read the lines data and attach them to the appropriate molecules in
    ! the catalog.
    text = ''
    nLines = 0
    do
      ! Find the next molecule heading
      do while ( index(text,'multicolumn') == 0 )
        read ( u, '(a)', end=8 ) text
      end do
      call getMolName
      ! Read the lines
      nLinesMol = 0
      do
        read ( u, '(a)', end=8 ) text
        if ( index(text,'multicolumn') /= 0 ) exit
        l = index(text,back)
        if ( l /= 0 ) then
          text(l:) = ''
        else
          l = len_trim(text) + 1
        end if
        if ( present(bands) ) then
          ! Find 15th &, which is before the list of bands
          i = 0
          do a = 1, 15
            i = index(text(i+1:l),'&') + i
          end do
          bands(nLines) = text(i+1:index(text(i+1:l),'&')+i-1)
        end if
        ! Replace TeX & separators by Fortran comma separators
        do i = 1, l
          if ( text(i:i) == '&' ) text(i:i) = ','
        end do
        nLines = nLines + 1       ! Total of all lines
        nLinesMol = nLinesMol + 1 ! Lines for this molecule
        catalog(mol)%lines(nLinesMol) = nLines
        read ( text(:l), * ) lines(nLines)%v0, lines(nLines)%el, &
          lines(nLines)%str, lines(nLines)%w, lines(nLines)%n, &
          lines(nLines)%ps, lines(nLines)%ns, lines(nLines)%delta, &
          lines(nLines)%n1, lines(nLines)%gamma, lines(nLines)%n2
        lines(nLines)%useYi = lines(nLines)%delta /= 0 .or. &
                              lines(nLines)%gamma /= 0
      end do
    end do
  8 continue
    close ( u )

    ! If there are no lines, allocate the lines component with zero size.
    ! Otherwise, when somebody looks at it, it will go bang.
    do mol = 1, size(catalog)
      if ( .not. associated(catalog(mol)%lines) ) then
        allocate ( catalog(mol)%lines(0), stat=stat )
        call test_allocate ( stat, moduleName, 'catalog(mol)%lines', (/ 1 /), &
          &  (/ 0 /), storage_size(catalog) / 8 )
      end if
    end do

    ! Destroy the cross-reference table
  9 continue
    deallocate ( TeX_name, xref_name )
    call test_deallocate ( stat, moduleName, 'TeX_name or XRef_name' )

  contains

    integer function FindMolecule ( TeX_text )
      character(*), intent(in) :: Tex_text
      integer :: i
      do i = 1, size(TeX_name)
        if ( index(trim(TeX_name(i)),trim(TeX_text)) /= 0 .and. &
           & index(trim(TeX_text),trim(TeX_name(i))) /= 0 ) then
          call getMoleculeIndex ( xref_name(i), findMolecule )
          return
        end if
      end do
      findMolecule = 0
    end function FindMolecule

    subroutine GetMolName
    ! Get the molecule name from the LaTeX table subheading, e.g.,
    ! \multicolumn{17}{c}{\dotfill $^{79}$BrO \dotfill} \\
      i = index(text,'dotfill')
      if ( i /= 0 ) then
        text = adjustl(text(i+7:))
        i = index(text,'dotfill')
        if ( i /= 0 ) text(i-1:) = ''
      end if
      mol = findMolecule ( text )
    end subroutine GetMolName

  end subroutine Read_Spectroscopy_TeX

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Read_Spectroscopy_TeX_m

! $Log$
! Revision 2.6  2015/03/28 02:05:30  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.5  2014/09/05 20:53:10  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.4  2014/08/05 01:16:43  vsnyder
! Use NEWUNIT= in OPEN statements instead of get_lun.  Eliminate dependence
! upon io_stuff.  Check for iostat == iostat_end instead of iostat /= 0.
! If there are no lines, allocate the lines component with zero size.
! Otherwise, when somebody looks at it, it will go bang.  Only calculate
! the base-10 logarithm of QLOG where it's > 0.
!
! Revision 2.3  2011/11/09 00:14:54  vsnyder
! Take base-10 logarithm of Qlog because that's what Create_Beta wants
!
! Revision 2.2  2011/11/01 22:11:51  vsnyder
! Remove the call to InitMolecules

! Revision 2.1  2011/11/01 20:50:06  vsnyder
! Initial commit
