! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2AUXData                 ! Data types for storing L2AUX data internally

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use Dump_0, only: DUMP
  use Hdf, only: DFACC_READ, DFNT_FLOAT64, SFCREATE, SFDIMID, SFSDSCALE, SFEND, &
    & SFENDACC, SFSTART, SFRDATA_F90, SFN2INDEX, SFSELECT, SFGINFO, &
    & SFGDINFO, SFSDMNAME, SFWDATA
  use intrinsic, only: LIT_INDICES, L_CHANNEL, L_GEODANGLE, L_LSBFREQUENCY, &
    & L_MAF, L_MIF, L_NONE, L_TIME, L_USBFREQUENCY
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
    & MLSMSG_ERROR, MLSMSG_WARNING
  use MLSSignals_m, only: GETMODULENAME, MODULES
  use MLSStrings, only: LINEARSEARCHSTRINGARRAY
  use Output_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
  use Tree, only: DUMP_TREE_NODE, SOURCE_REF

  implicit none

  ! Externals

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2AUX data.
  ! L2AUX data is a rank three array, with the first index being a channel like
  ! index, the second a height like index, and the last an instance index.  It
  ! is directly conformable to a vector quantity.
  ! If a dimension is absent, it is flagged appropriately but the array is still
  ! rank 3.  However, when an HDF file is written or read missing dimensions
  ! are swallowed.

  ! For Announce_Error
    integer :: ERROR

  integer, parameter :: L2AUXRANK=3     ! Dimensionality of L2AUXData_T%values

  ! This datatype describes a dimension for an L2AUX quantity

  type L2AUX_Dimension_T
     integer :: NOVALUES        ! Length of this dimension
     integer :: DIMENSIONFAMILY ! What is this dimension
     real(r8), dimension(:), pointer :: VALUES=>NULL() ! (noValues)
  end type L2AUX_Dimension_T

  ! This datatype describes an l2aux quantity itself.
  ! The dimensions will typically be ordered as follows:
  ! [Channel or frequency], MIF, [MAF or time or geodAngle]

  type L2AUXData_T
    integer :: NAME                     ! String index of name to be output
    integer :: INSTRUMENTMODULE         ! From source vector
    logical :: MINORFRAME               ! Is this a minor frame quantity
    logical :: MAJORFRAME               ! Is this a major frame quantity
    ! The dimensions for the quantity
    type (L2AUX_Dimension_T), dimension(L2AUXRank) :: DIMENSIONS
    ! The values of the quantitiy
    real(r8), pointer, dimension(:,:,:) :: VALUES=>NULL()
  end type L2AUXData_T

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------  SetupNewL2AUXRecord   -----
  subroutine SetupNewL2AUXRecord ( dimensionFamilies, dimSizes, dimStarts, l2aux )

    ! This first routine sets up the arrays for an l2aux datatype.
    ! The user supplies a set of three dimensionFamilies (e.g. l_maf)
    ! Quantities can have upto three valid dimensions.  l_none can be used
    ! to indicate later dimensions are invalid.

    ! Dummy arguments
    integer, dimension(L2AUXRank), intent(in) :: dimensionFamilies
    integer, dimension(L2AUXRank), intent(in) :: dimSizes
    integer, dimension(L2AUXRank), intent(in) :: dimStarts
    type (L2AUXData_T), intent(out) :: l2aux

    ! Local variables
    integer :: dimIndex
    integer :: status
    integer, dimension(L2AUXRank) :: dimEnds

    ! Fill the dimensions data structure
    l2aux%dimensions%dimensionFamily = dimensionFamilies
    l2aux%dimensions%noValues = dimSizes

    dimEnds = dimStarts + max(1,dimSizes) - 1

    ! Allocate the values for each dimension
    do dimIndex = 1, L2AUXRank
      if ( dimensionFamilies(dimIndex)/=L_None ) then
        allocate (l2aux%dimensions(dimIndex)%values( &
          & dimStarts(dimIndex):dimEnds(dimIndex)), &
          & STAT=status)
        if ( status/=0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate // "l2aux dimension values" )
      else
        l2aux%dimensions(dimIndex)%noValues=1
      end if
    end do

    ! Allocate the values for the data itself

    allocate ( l2aux%values( &
      & dimStarts(1):dimEnds(1), &
      & dimStarts(2):dimEnds(2), &
      & dimStarts(3):dimEnds(3)), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate// "l2aux values" )

  end subroutine SetupNewL2AUXRecord
    
  !----------------------------------------  DestroyL2AUXContents  -----
  subroutine DestroyL2AUXContents ( l2aux )

  ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: l2aux
    ! Local variables
    integer :: dim
    ! Executable code
    
    do dim=1,L2AUXRank
      if (l2aux%dimensions(dim)%dimensionFamily /= L_None) &
        call deallocate_test(l2aux%dimensions(dim)%values, &
        "l2aux%dimensions", ModuleName)
    end do
    call deallocate_test(l2aux%values,"l2aux%values",ModuleName ) 
  end subroutine DestroyL2AUXContents

  !--------------------------------------  ExpandL2AUXDataInPlace  -----
  subroutine ExpandL2AUXDataInPlace ( l2aux, newSize )

  ! This subroutine expands an L2AUXData_T in place, allowing the user to
  ! add more `profiles' to it.  Note that the `profile' dimension is the last
  ! one.

    ! Dummy arguments
    type (L2AUXData_T), intent(inout) :: L2aux
    integer, intent(in) :: NewSize

    ! Local variables
    integer :: OldSize
    ! The following are temporary arrays for copying data around
    real (r8), dimension(:), pointer :: Temp1D
    real (r8), dimension(:,:,:), pointer :: Temp3D

    ! Executable code

    nullify ( temp1D, temp3D )

    if (l2aux%dimensions(3)%dimensionFamily==L_None) &
      call MLSMessage (MLSMSG_Error, ModuleName, &
        & "This l2aux is not expandable")

    ! Now see how long this is
    oldSize = l2aux%dimensions(3)%noValues
    ! Do a sanity check
    if ( newSize<oldSize ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "This l2aux is getting smaller not bigger" )

    ! Now expand this dimension
    temp1D => l2aux%dimensions(3)%values
    ! Nullify old one so allocate_test doesn't clobber temp1D
    nullify ( l2aux%dimensions(3)%values )
    call allocate_test(l2aux%dimensions(3)%values,newSize, &
      & 'New l2aux%dimensions(3)%values', ModuleName)
    l2aux%dimensions(3)%noValues = newSize
    l2aux%dimensions(3)%values(1:oldSize)=temp1D
    ! Now we can loose the old values
    call deallocate_test ( temp1D, "temp1D", ModuleName )

    ! Now expand the data in this dimension, save old field
    temp3d => l2aux%values
    ! Need to nullify old one, else next call will clobber it
    ! loosing us temp3d
    nullify ( l2aux%values )
    call allocate_test ( l2aux%values, &
      & max(1, l2aux%dimensions(1)%noValues), &
      & max(1, l2aux%dimensions(2)%noValues), &
      & newSize, "l2aux%values", ModuleName )
    l2aux%values(:,:,1:oldSize) = temp3d

    ! Now we can loose temp3d
    call deallocate_test ( temp3d, "temp3d", ModuleName )

  end subroutine ExpandL2AUXDataInPlace

  !------------------------------------------  AddL2AUXToDatabase  -----
  integer function AddL2AUXToDatabase ( DATABASE, ITEM )

  ! This subroutine adds an l2aux data type to a database of said types,
  ! creating a new database if it doesn't exist.  The result value is
  ! the length of the database -- where L2aux is put.

    ! Dummy arguments
    type (L2AUXData_T), dimension(:), pointer :: DATABASE
    type (L2AUXData_T), intent(in) :: ITEM

    ! local variables
    type (L2AUXData_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddL2AUXToDatabase = newSize
  end function AddL2AUXToDatabase

  ! ---------------------------------------  DestroyL2AUXDatabase  -----
  subroutine DestroyL2AUXDatabase ( DATABASE )

  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (L2AUXData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2auxIndex, status

    if ( associated(database) ) then
      do l2auxIndex = 1, size(database)
        call DestroyL2AUXContents ( database(l2auxIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Warning,ModuleName, &
        & MLSMSG_DeAllocate // "database")
    end if
  end subroutine DestroyL2AUXDatabase

  ! ------------------------------------------ Dump_L2AUX ------------

  subroutine Dump_L2AUX ( L2aux, Name )

    ! Dummy arguments
    type (l2auxData_T), intent(in) ::          L2AUX(:)
    character(len=*), intent(in), optional :: Name

    ! Local variables
    integer :: i, dim
    
    if ( present(name) ) call output ( name, advance='yes' )
    do i = 1, size(l2aux)
      call output ( 'L2AUX Data: ')
      call display_string ( l2aux(i)%name )
      call output ( '    instrumentmodule: ')
!      call display_string ( l2aux(i)%instrumentmodule )
      call display_string ( modules(l2aux(i)%instrumentmodule)%name, advance='yes' ) 
      call output ( '    (its index): ')
      call output ( l2aux(i)%instrumentmodule, advance='no')
      call output ( ' ', advance='yes')
      call output ( '  Minor Frame? (t/f): ')
      call output ( l2aux(i)%minorframe, advance='no')
      call output ( '  Major Frame? (t/f): ')
      call output ( l2aux(i)%majorframe, advance='yes')
      do dim=1, l2auxrank
        call output ( '  dimension: ')
        call output ( dim )
        call output ( '           ')
        if ( associated(l2aux(i)%dimensions(dim)%values) ) then
          call output ( '  nValues: ')
          call output ( l2aux(i)%dimensions(dim)%novalues, 3, advance='no')
          call output ( '           ')
          call output ( '  dimension family: ')
          call output ( l2aux(i)%dimensions(dim)%dimensionfamily, 3, advance='yes')
          call dump ( l2aux(i)%dimensions(dim)%values, 'dim values:' )
         else
        call output ( ' is not associated', advance='yes')
         endif
      enddo
      call dump ( l2aux(i)%values, 'values:' )
 
    end do
  end subroutine Dump_L2AUX
    
  !------------------------------------------------ ReadL2AUXData ------------
  subroutine ReadL2AUXData(sd_id, quantityname, l2aux, firstProf, lastProf)
    
    ! This routine reads an l2aux file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    character (LEN=*), intent(IN) :: quantityname ! Name of L2AUX quantity = sdname in writing routine
    integer, intent(IN) :: sd_id ! Returned by sfstart before calling us
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2AUXData_T ), intent(OUT) :: l2aux ! Result

    ! Parameters

    character (LEN=*), parameter :: SZ_ERR = 'Failed to get size of &
         &dimension '
    character (LEN=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (LEN=*), parameter :: MLSMSG_l2auxRead = 'Unable to read l2aux &
                                                     &field:'
    integer, parameter :: MAXRANK = 3
    integer, parameter :: MAXDIMSIZES = 300
    logical, parameter :: CHECKDIMSIZES = .true.	! .TRUE. only while debugging

    ! Functions

    !INTEGER, EXTERNAL :: swattach, swdetach, swdiminfo, swinqdims, swrdfld

    ! Variables

    character (LEN=480) :: msr

    integer :: sds_index, sds_id, rank, data_type, num_attrs, dim, dim_id
    integer :: dim_sizes(MAXRANK), dim_size1
    integer :: dim_families(MAXRANK)
    character (LEN=len(quantityname)) :: sds_name
    character (LEN=132) :: dim_name
    character (LEN=1)                  :: dim_char

    integer :: nDims, status
    integer :: start(3), stride(3), edge(3), dims(3)

    logical :: firstCheck, lastCheck

    ! Attach to the file for reading

    ! find SD data set identifier
    sds_index = sfn2index(sd_id, quantityname)
    if (sds_index == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get sds_index.')
         
    sds_id = sfselect(sd_id, sds_index)
    if (sds_id == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get sds_id.')
         
    status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, &
    & num_attrs)

    if (status == -1) then
       call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         & get sf info.')
    elseif (sds_name /= quantityname) then
       call MLSMessage(MLSMSG_Error, ModuleName, 'quantityname &
         & fails to match sf info.')
    endif

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    ! Uncertain what to do with those just yet
    ! Now find dimension family of dimension; e.g., MAF
    dim_size1 = 0
     do dim=1, rank

    	write(dim_char, '(I1)') dim
    	dim_id = sfdimid(sds_id, dim-1)		! dim starts from 0
        if(dim_id == -1) then

           msr = 'Failed to &
           & get dim_id for dim index number ' // dim_char
           call MLSMessage(MLSMSG_Error, ModuleName, msr)
        else
            status = sfgdinfo(dim_id, dim_name, dim_size1, data_type, &
            & num_attrs)

            if(status == -1) then
                  msr = 'Failed to &
                  & get dim_info for dim index number ' // dim_char
                  call MLSMessage(MLSMSG_Error, ModuleName, msr)
            else
                dim_families(dim) = dim_size1
                if(dim_families(dim) == 0) then
                     msr = 'Failed to &
                     & find ' //dim_name // ' among L2AuxDimNames'
                     call MLSMessage(MLSMSG_Error, ModuleName, msr)
                endif
            endif
         endif
    enddo
    ! Allocate result

    call SetupNewl2auxRecord ( dim_families, dim_sizes, (/1,1,1/), l2aux )

    ! Allocate temporary arrays

    ! Read the SD
    start = 0
    stride = 1
    status = sfrdata_f90(sds_id, start, stride, dim_sizes, l2aux%values)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         & write SD.')

    ! Deallocate local variables


    ! Terminate access to the data set

    status = sfendacc(sds_id)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &end access to sds_id after reading.')

    !  After reading, detach from hdf interface

    status = sfend(sd_id)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from SD file after reading.')

  end subroutine ReadL2AUXData

  !----------------------------------------------------- WriteL2AUXData ------

  subroutine WriteL2AUXData(l2aux, l2FileHandle, sdName)
    type (L2AUXData_T), intent(in) :: L2AUX
    integer, intent(in) :: L2FILEHANDLE
    character (len=*), optional, intent(in) :: SDNAME ! Defaults to l2aux%name

    ! Local variables
    integer :: NODIMENSIONSUSED         ! No. real dimensions
    integer, dimension(:), pointer :: DIMSIZES ! Size of each dimension
    integer :: sdId                   ! sd id from sfcreate
    logical, dimension(L2AUXRank) :: GOODDIM ! Not L_None
    integer :: dimensionInData     ! Index
    integer :: dimensionInFile          ! Index
    integer :: dimID                    ! ID for HDF
    character (len=132) :: NAMESTRING   ! Name for sd
    character (len=132) :: DIMNAME      ! Name for dimension
    integer :: status                   ! Flag
    integer, parameter, dimension(L2AUXRank) :: stride = (/ 1, 1, 1/)
    integer, parameter, dimension(L2AUXRank) :: start = (/ 0, 0, 0/)


    ! Executable code
    nullify ( dimSizes )
    error = 0

    if (present(sdName)) then
      nameString=sdName
    else
      call get_string ( l2aux%name, nameString, strip=.true. )
    endif

    goodDim=l2aux%dimensions%dimensionFamily /= L_None
    noDimensionsUsed=COUNT(goodDim)
    call allocate_test(dimSizes,noDimensionsUsed,'dimSizes',ModuleName)
    dimSizes=PACK(l2aux%dimensions%noValues, goodDim)

    ! Create the sd within the file
    sdId= SFcreate ( l2FileHandle, nameString, DFNT_FLOAT64, &
      & noDimensionsUsed, dimSizes)

    ! Now define the dimensions
    dimensionInFile=0
    do dimensionInData=1,L2AUXRank
      if (l2aux%dimensions(dimensionInData)%dimensionFamily &
        &    /= L_None) then
        dimID=SFDIMID(sdId,dimensionInFile)
        ! Construct dimension name. For minor frame quantities MIF and MAF
        ! names are global, otherwise specific
        if ( (dimensionInData > 1) .and. (l2aux%minorFrame) ) then
          call GetModuleName( l2aux%instrumentModule, dimName)
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
        else
          call get_string (l2aux%name, dimName, strip=.true. )
          if (len_trim(dimName) < len(dimName)) dimName=TRIM(dimName)//'.'
        endif
        call get_string (lit_indices(l2aux%dimensions(dimensionInData)%dimensionFamily), &
          & dimName(LEN_TRIM(dimName)+1:))
        ! Write dimension name
        status=SFSDMName(dimID,TRIM(dimName))
        if ( status /= 0 ) then
		  		call output("dim name: ")
		  		call output(TRIM(dimName), advance='yes')
		  		call announce_error (  0, &
          & "Error setting dimension name to SDS l2aux file:")
!		  		call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error setting dimension name to SDS l2aux file:")
	endif
        ! Write dimension scale
        status=SFSDScale(dimID, dimSizes(dimensionInFile+1), DFNT_FLOAT64,&
          & l2aux%dimensions(dimensionInData)%values)
        if ( status /= 0 ) then
		  		call output("dimID: ")
		  		call output(dimID, advance='yes')
		      call announce_error ( 0, &
          & "Error writing dimension scale in l2auxFile:" )
!		      call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error writing dimension scale in l2auxFile:" )
	endif
        dimensionInFile=dimensionInFile+1
      endif
    end do

    ! Now write the data
    status= SFWData(sdId, start(1:noDimensionsUsed), &
      & stride(1:noDimensionsUsed), dimSizes, l2aux%values)
    if ( status /= 0 ) then
	   call announce_error (0,&
      & "Error writing SDS data to  l2aux file:  " )
!	   call MLSMessage ( MLSMSG_Error, ModuleName,&
!      & "Error writing SDS data to  l2aux file:  " )
    endif

    call Deallocate_Test(dimSizes,"dimSizes",ModuleName)
    
    ! Terminate access to sd
    status = sfendacc(sdId)

  end subroutine WriteL2AUXData

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, full_message, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
	character(LEN=*), intent(in)    :: full_message
    integer, intent(in), optional :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
!    call print_source ( source_ref(where) )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ' OutputAndClose complained: ' )


    call output ( " Caused the following error: ", advance='yes', &
      & from_where=ModuleName )
    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      select case ( code )
      end select
    end if
  end subroutine ANNOUNCE_ERROR

!=============================================================================
end module L2AUXData
!=============================================================================

!
! $Log$
! Revision 2.20  2001/10/08 23:41:27  pwagner
! Improved dump routines
!
! Revision 2.19  2001/10/05 23:32:27  pwagner
! Added majorframe to data type; trimmed unused stuff
!
! Revision 2.18  2001/08/06 18:35:24  pwagner
! Added dump_l2aux
!
! Revision 2.17  2001/07/11 20:50:46  dwu
! fix problem in readl2auxdata
!
! Revision 2.16  2001/05/30 23:53:31  livesey
! Changed for new version of L1BData
!
! Revision 2.15  2001/05/12 00:18:40  livesey
! Tidied up issues with array bounds etc.
!
! Revision 2.14  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.13  2001/05/02 22:24:20  pwagner
! Removed SDPToolkit use
!
! Revision 2.12  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.11  2001/04/12 22:19:33  vsnyder
! Improved an error message
!
! Revision 2.10  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.9  2001/04/07 00:14:27  pwagner
! Added announce_error
!
! Revision 2.8  2001/03/15 18:42:29  livesey
! Removed quotes from dimension name prefixes
!
! Revision 2.7  2001/03/08 02:20:12  livesey
! Added strip argument to a call to get_string
!
! Revision 2.6  2001/03/06 22:40:47  livesey
! Working version
!
! Revision 2.5  2001/02/14 23:41:33  livesey
! Removed irrelevant numProfs argument
!
! Revision 2.4  2001/01/03 00:46:19  pwagner
! Changed sfgetinfo to sfginfo
!
! Revision 2.3  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
! Revision 2.2  2000/12/04 21:48:29  pwagner
! ReadL2AUXData completed
!
! Revision 2.1  2000/12/02 01:12:00  pwagner
! Added ReadL2AUXData
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

