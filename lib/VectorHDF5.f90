! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module VectorHDF5
  ! This module reads and writes vectors (complete with quantity templates etc.
  ! to and from HDF5 files.

  implicit none
  private

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  public :: WriteVectorAsHDF5, ReadVectorFromHDF5

contains ! ========================================= Module procedures =======

  subroutine WriteQuantityTemplateAsHDF5 ( location, qt, index )
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use QuantityTemplates, only: QUANTITYTEMPLATE_T
    use String_table, only: GET_STRING
    use MLSHDF5, only: MAKEHDF5ATTRIBUTE, WRITESTRINGINDEXASHDF5ATTRIBUTE, &
      & WRITELITINDEXASHDF5ATTRIBUTE, SAVEASHDF5DS
    use HDF5, only: H5GCREATE_F, H5GCLOSE_F

    ! Arguments
    integer, intent(in) :: LOCATION     ! HDF location
    type(QuantityTemplate_T), intent(in) :: QT
    integer, intent(in) :: INDEX        ! Index into parent vector

    ! Local variables
    character(len=132) :: QNAME         ! Name of quantity
    integer :: QID                      ! HDF5 ID of quantity
    integer :: STATUS                   ! Flag from HDF5
    ! Executable code
    ! Create a group for the quantity
    call get_string ( qt%name, qName )
    call h5gCreate_f ( location, trim(qName), qID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for quantity ' // trim(qName) )

    ! Write the low level stuff as attributes
    call WriteStringIndexAsHDF5Attribute ( qID, 'name',qt%name )
    call MakeHDF5Attribute ( qID, 'id', qt%id )
    call MakeHDF5Attribute ( qID, 'index', index )
    call WriteLitIndexAsHDF5Attribute ( qID, 'quantityType', qt%quantityType )
    call MakeHDF5Attribute ( qID, 'noChans', qt%noChans )
    call MakeHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
    call MakeHDF5Attribute ( qID, 'noInstances', qt%noInstances )
    call MakeHDF5Attribute ( qID, 'instanceLen', qt%instanceLen )
    call MakeHDF5Attribute ( qID, 'coherent', qt%coherent )
    call MakeHDF5Attribute ( qID, 'stacked', qt%stacked )
    call MakeHDF5Attribute ( qID, 'regular', qt%regular )
    call MakeHDF5Attribute ( qID, 'minorFrame', qt%minorFrame )
    call MakeHDF5Attribute ( qID, 'majorFrame', qt%majorFrame )
    call MakeHDF5Attribute ( qID, 'logBasis', qt%logBasis )
    call MakeHDF5Attribute ( qID, 'minValue', qt%minValue )
    call MakeHDF5Attribute ( qID, 'noInstancesLowerOverlap', qt%noInstancesLowerOverlap )
    call MakeHDF5Attribute ( qID, 'noInstancesUpperOverlap', qt%noInstancesUpperOverlap )
    call WriteLitIndexAsHDF5Attribute ( qID, 'verticalCoordinate', qt%verticalCoordinate )
    call MakeHDF5Attribute ( qID, 'badValue', qt%badValue )
    call WriteLitIndexAsHDF5Attribute ( qID, 'unit', qt%unit )
    call MakeHDF5Attribute ( qID, 'scaleFactor', qt%scaleFactor )

    call WriteLitIndexAsHDF5Attribute ( qID, 'frequencyCoordinate', qt%frequencyCoordinate )
    call MakeHDF5Attribute ( qID, 'lo', qt%lo )
    call MakeHDF5Attribute ( qID, 'signal', qt%signal )
    call MakeHDF5Attribute ( qID, 'sideband', qt%sideband )
    call MakeHDF5Attribute ( qID, 'instrumentModule', qt%instrumentModule )
    call MakeHDF5Attribute ( qID, 'radiometer', qt%radiometer )
    call MakeHDF5Attribute ( qID, 'molecule', qt%molecule )

    ! Now write the coordinates
    call SaveAsHDF5DS ( qID, 'surfs', qt%surfs )
    call SaveAsHDF5DS ( qID, 'phi', qt%phi )
    call SaveAsHDF5DS ( qID, 'geodLat', qt%geodLat )
    call SaveAsHDF5DS ( qID, 'lon', qt%lon )
    call SaveAsHDF5DS ( qID, 'time', qt%time )
    call SaveAsHDF5DS ( qID, 'solarTime', qt%solarTime )
    call SaveAsHDF5DS ( qID, 'solarZenith', qt%solarZenith )
    call SaveAsHDF5DS ( qID, 'losAngle', qt%losAngle )

    ! These may or may not be present
    if ( associated ( qt%mafIndex ) ) &
      & call SaveAsHDF5DS ( qID, 'mafIndex', qt%mafIndex )
    if ( associated ( qt%mafCounter ) ) &
      & call SaveAsHDF5DS ( qID, 'mafCounter', qt%mafCounter )
    if ( associated ( qt%frequencies ) ) &
      & call SaveAsHDF5DS ( qID, 'frequencies', qt%frequencies )
    if ( associated ( qt%surfIndex ) ) &
      & call SaveAsHDF5DS ( qID, 'surfIndex', qt%surfIndex )
    if ( associated ( qt%chanIndex ) ) &
      & call SaveAsHDF5DS ( qID, 'chanIndex', qt%chanIndex )

    ! Close the group
    call h5gClose_f ( qID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close quantity group for ' // trim(qName) )
  end subroutine WriteQuantityTemplateAsHDF5

  ! ----------------------------------------- WriteVectorAsHDF5 -----
  subroutine WriteVectorAsHDF5 ( location, name, vector )
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use String_table, only: GET_STRING
    use VectorsModule, only: VECTOR_T
    use MLSHDF5, only: MAKEHDF5ATTRIBUTE, WRITESTRINGINDEXASHDF5ATTRIBUTE, &
      & SAVEASHDF5DS
    use HDF5, only: H5GCREATE_F, H5GCLOSE_F
    ! Arguments
    integer, intent(in) :: LOCATION     ! HDF location
    character(len=*), intent(in) :: NAME ! Name of group to put it in
    type(Vector_T), intent(in) :: VECTOR

    ! Local variables
    integer :: VID                      ! HDF5 ID of vector
    integer :: STATUS                   ! Flag from HDF5
    integer :: GID                      ! HDF5 ID of subgroup
    integer :: Q                        ! Loop counter
    character(len=132) :: WORD          ! A name for something

    ! Executable code
    ! Create group for vector
    call h5gCreate_f ( location, trim(name), vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group for vector ' // trim(name) )

    ! Annotate it
    call MakeHDF5Attribute ( vID, 'noQuantities', size(vector%quantities) )
    call MakeHDF5Attribute ( vID, 'id', vector%template%id )
    call MakeHDF5Attribute ( vID, 'totalInstances', vector%template%totalInstances )
    call MakeHDF5Attribute ( vID, 'totalElements', vector%template%totalElements )
    call WriteStringIndexAsHDF5Attribute ( vID, 'templateName', vector%template%name )

    ! Stick quantity templates in one group
    call h5gCreate_f( vId, 'templates', gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create templates group for vector ' // trim(name) )
    ! Write out the quantity templates
    do q = 1, size ( vector%quantities )
      call WriteQuantityTemplateAsHDF5 ( gid, vector%quantities(q)%template, q )
    end do
    ! Close the group
    call h5gClose_f ( gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close templates group for vector ' // trim(name) )

    ! Create a group for the quantity values
    call h5gCreate_f( vId, 'values', gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create values group for vector ' // trim(name) )
    ! Write out the quantity values
    do q = 1, size ( vector%quantities )
      call get_string ( vector%quantities(q)%template%name, word, &
        & strip=.true., noError=.true. )
      call SaveAsHDF5DS ( gid, trim(word), vector%quantities(q)%values )
    end do
    ! Close the group
    call h5gClose_f ( gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close values group for vector ' // trim(name) )

    ! Create a group for the quantity masks
    call h5gCreate_f( vId, 'masks', gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create masks group for vector ' // trim(name) )
    ! Write out the quantity values
    do q = 1, size ( vector%quantities )
      if ( associated ( vector%quantities(q)%mask ) ) then
        call get_string ( vector%quantities(q)%template%name, word, &
          & strip=.true., noError=.true. )
        call SaveAsHDF5DS ( gid, trim(word), vector%quantities(q)%mask )
      end if
    end do
    ! Close the group
    call h5gClose_f ( gid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close masks group for vector ' // trim(name) )

    ! Finished, close the group
    call h5gClose_f ( vid, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group for vector ' // trim(name) )

  end subroutine WriteVectorAsHDF5

  ! ----------------------------------- ReadQuantityTemplateFromHDF5 ---
  subroutine ReadQuantityTemplateFromHDF5 ( location, name, qt )
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use HDF5, only: H5GOPEN_F, H5GCLOSE_F
    use QuantityTemplates, only: QUANTITYTEMPLATE_T, SETUPNEWQUANTITYTEMPLATE
    use MLSHDF5, only: LOADFROMHDF5DS, GETHDF5ATTRIBUTE, &
      & READSTRINGINDEXFromHDF5Attr, READLITINDEXFromHDF5Attr
    ! Arguments
    integer, intent(in) :: LOCATION     ! HDF5 location of parent group
    character(len=*), intent(in) :: NAME ! Name of quantity template
    type(QuantityTemplate_T), intent(out) :: QT ! Template for quantity
    
    ! Local variables
    integer :: NOINSTANCES              ! Dimension
    integer :: NOSURFS                  ! Dimension
    integer :: NOCHANS                  ! Dimension
    integer :: INSTANCELEN              ! Dimension
    logical :: COHERENT                 ! Flag
    logical :: STACKED                  ! Flag
    logical :: REGULAR                  ! Flag
    logical :: MAJORFRAME               ! Flag
    logical :: MINORFRAME               ! Flag
    integer :: QID                     ! HDF5 ID of quantity template
    integer :: STATUS                   ! Flag from HDF

    ! Executable code
    ! Attach to the group
    call h5gOpen_f ( location, trim(name), qID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group for quantity template '//trim(name) )

    ! Get the fundamentals
    call GetHDF5Attribute ( qID, 'noChans', noChans )
    call GetHDF5Attribute ( qID, 'noSurfs', noSurfs )
    call GetHDF5Attribute ( qID, 'noInstances', noInstances )
    call GetHDF5Attribute ( qID, 'instanceLen', instanceLen )
    call GetHDF5Attribute ( qID, 'coherent', coherent )
    call GetHDF5Attribute ( qID, 'stacked', stacked )
    call GetHDF5Attribute ( qID, 'regular', regular )
    call GetHDF5Attribute ( qID, 'minorFrame', minorFrame )
    call GetHDF5Attribute ( qID, 'majorFrame', majorFrame )
    
    call SetupNewQuantityTemplate ( qt, noChans=noChans, noSurfs=noSurfs, &
      & noInstances=noInstances, coherent=coherent, stacked=stacked, &
      & regular=regular, minorFrame=minorFrame, majorFrame=majorFrame, &
      & instanceLen=instanceLen )

    ! Get the remaining stuff
    ! Write the low level stuff as attributes
    call ReadStringIndexFromHDF5Attr ( qID, 'name', qt%name )
    call GetHDF5Attribute ( qID, 'id', qt%id )
    call ReadLitIndexFromHDF5Attr ( qID, 'quantityType', qt%quantityType )
    call GetHDF5Attribute ( qID, 'noChans', qt%noChans )
    call GetHDF5Attribute ( qID, 'noSurfs', qt%noSurfs )
    call GetHDF5Attribute ( qID, 'noInstances', qt%noInstances )
    call GetHDF5Attribute ( qID, 'instanceLen', qt%instanceLen )
    call GetHDF5Attribute ( qID, 'coherent', qt%coherent )
    call GetHDF5Attribute ( qID, 'stacked', qt%stacked )
    call GetHDF5Attribute ( qID, 'regular', qt%regular )
    call GetHDF5Attribute ( qID, 'minorFrame', qt%minorFrame )
    call GetHDF5Attribute ( qID, 'majorFrame', qt%majorFrame )
    call GetHDF5Attribute ( qID, 'logBasis', qt%logBasis )
    call GetHDF5Attribute ( qID, 'minValue', qt%minValue )
    call GetHDF5Attribute ( qID, 'noInstancesLowerOverlap', qt%noInstancesLowerOverlap )
    call GetHDF5Attribute ( qID, 'noInstancesUpperOverlap', qt%noInstancesUpperOverlap )
    call ReadLitIndexFromHDF5Attr ( qID, 'verticalCoordinate', qt%verticalCoordinate )
    call GetHDF5Attribute ( qID, 'badValue', qt%badValue )
    call ReadLitIndexFromHDF5Attr ( qID, 'unit', qt%unit )
    call GetHDF5Attribute ( qID, 'scaleFactor', qt%scaleFactor )
    call ReadLitIndexFromHDF5Attr ( qID, 'frequencyCoordinate', qt%frequencyCoordinate )
    call GetHDF5Attribute ( qID, 'lo', qt%lo )
    call GetHDF5Attribute ( qID, 'signal', qt%signal )
    call GetHDF5Attribute ( qID, 'sideband', qt%sideband )
    call GetHDF5Attribute ( qID, 'instrumentModule', qt%instrumentModule )
    call GetHDF5Attribute ( qID, 'radiometer', qt%radiometer )
    call GetHDF5Attribute ( qID, 'molecule', qt%molecule )

    ! Now write the coordinates
    call LoadFromHDF5DS ( qID, 'surfs', qt%surfs )
    call LoadFromHDF5DS ( qID, 'phi', qt%phi )
    call LoadFromHDF5DS ( qID, 'geodLat', qt%geodLat )
    call LoadFromHDF5DS ( qID, 'lon', qt%lon )
    call LoadFromHDF5DS ( qID, 'time', qt%time )
    call LoadFromHDF5DS ( qID, 'solarTime', qt%solarTime )
    call LoadFromHDF5DS ( qID, 'solarZenith', qt%solarZenith )
    call LoadFromHDF5DS ( qID, 'losAngle', qt%losAngle )

    ! These may or may not be present
    if ( associated ( qt%mafIndex ) ) &
      & call LoadFromHDF5DS ( qID, 'mafIndex', qt%mafIndex )
    if ( associated ( qt%mafCounter ) ) &
      & call LoadFromHDF5DS ( qID, 'mafCounter', qt%mafCounter )
    if ( associated ( qt%frequencies ) ) &
      & call LoadFromHDF5DS ( qID, 'frequencies', qt%frequencies )
    if ( associated ( qt%surfIndex ) ) &
      & call LoadFromHDF5DS ( qID, 'surfIndex', qt%surfIndex )
    if ( associated ( qt%chanIndex ) ) &
      & call LoadFromHDF5DS ( qID, 'chanIndex', qt%chanIndex )

    ! Close the group
    call h5gClose_f ( qID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close quantity group for ' // trim(name) )
  end subroutine ReadQuantityTemplateFromHDF5

  ! -------------------------------------------- ReadVectorFromHDF5
  subroutine ReadVectorFromHDF5 ( location, name, vector, &
    & quantities, noQuantityTemplates )
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use VectorsModule, only: VECTOR_T, VECTORTEMPLATE_T, CREATEVECTOR, CREATEMASK
    use MoreTree, only: GETSTRINGINDEXFROMSTRING
    use String_Table, only: GET_STRING
    use QuantityTemplates, only: NULLIFYQUANTITYTEMPLATE, QUANTITYTEMPLATE_T, &
      & ADDQUANTITYTEMPLATETODATABASE, INFLATEQUANTITYTEMPLATEDATABASE
    use MLSHDF5, only: ISHDF5DSPRESENT, GETHDF5ATTRIBUTE, LOADFROMHDF5DS, &
      & READSTRINGINDEXFromHDF5Attr
    use HDF5, only: H5GOPEN_F, H5GCLOSE_F, H5GGET_OBJ_INFO_IDX_F

    ! Arguments
    integer, intent(in) :: LOCATION     ! HDF5 location
    character(len=*), intent(in) :: NAME ! Name of group vector is in
    type (Vector_T), intent(out) :: VECTOR ! Resulting vector
    type (QuantityTemplate_T), dimension(:), pointer :: QUANTITIES
    integer, intent(inout), optional :: NOQUANTITYTEMPLATES

    ! Local parameters 
    integer, parameter :: DATABASEINFLATION = 500

    ! Local variables
    integer :: VID                      ! Group ID for vector
    integer :: GID                      ! A subgroup id
    integer :: QID                      ! A quantity template id
    integer :: STATUS                   ! Flag from HDF
    integer :: NOQUANTITIES             ! Number of vector quantities
    integer :: OBJTYPE                  ! From HDF5
    integer :: MYNOQUANTITYTEMPLATES    ! Number of quantity templates
    integer :: Q                        ! Loop counter
    integer :: I                        ! Index of quantity
    character (len=64), pointer, dimension(:) :: QUANTITYNAMES ! Names of quantities
    character (len=64) :: WORD          ! Name of quantity

    type (VectorTemplate_T) :: VT       ! Vector template
    type (QuantityTemplate_T) :: QT     ! QuantityTemplate
    
    ! Executable code
    if ( associated ( quantities ) ) then
      myNoQuantityTemplates = size ( quantities )
      if ( present ( noQuantityTemplates ) ) &
        & myNoQuantityTemplates = noQuantityTemplates
    else
      myNoQuantityTemplates = 0
    end if

    ! Attach to group for vector
    call h5gOpen_f ( location, trim(name), vID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group for vector ' // trim(name) )

    ! Get the information we need.
    call GetHDF5Attribute ( vID, 'noQuantities', vt%noQuantities )
    call GetHDF5Attribute ( vID, 'id', vt%id )
    call GetHDF5Attribute ( vID, 'totalInstances', vt%totalInstances )
    call GetHDF5Attribute ( vID, 'totalElements', vt%totalElements )
    call ReadStringIndexFromHDF5Attr ( vID, 'templateName', vt%name )

    ! Open the quantity templates
    call h5gOpen_f ( vId, 'templates', gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open templates group for vector ' // trim(name) )

    ! Work out the quantity names in order
    nullify ( quantityNames )
    call Allocate_test ( quantityNames, vt%noQuantities, 'quantityNames', ModuleName )
    call Allocate_test ( vt%quantities, vt%noQuantities, 'templates', ModuleName )

    do q = 1, vt%noQuantities
      call h5gget_obj_info_idx_f ( vId, 'templates', q-1, word, &
        & objType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to identify a quantity within '//trim(name) )
      call h5gOpen_f ( gId, trim(word), qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open quantity '//trim(word)//' in vector '//trim(name) )
      call GetHDF5Attribute ( qId, 'index', i )
      quantityNames ( i ) = word
      call h5gClose_f ( qId, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close quantity '//trim(word)//' in vector '//trim(name) )
    end do

    ! Now read the quantity templates
    do q = 1, vt%noQuantities
      ! First make sure we don't tread on previous quantities
      ! I know the read routine is intent(out) but you never know
      call NullifyQuantityTemplate ( qt )
      call ReadQuantityTemplateFromHDF5 ( gId, trim(quantityNames(q)), qt )
      ! Now add this to our database
      if ( present ( noQuantityTemplates ) ) then
        ! Do this by inflation
        if ( myNoQuantityTemplates > size ( quantities ) ) &
          & myNoQuantityTemplates = InflateQuantityTemplateDatabase ( &
          & quantities, DatabaseInflation )
        myNoQuantityTemplates = myNoQuantityTemplates + 1
        quantities ( myNoQuantityTemplates ) = qt
      else
        myNoQuantityTemplates = AddQuantityTemplateToDatabase ( &
          & quantities, qt )
      end if
      vt%quantities ( q ) = myNoQuantityTemplates
    end do

    ! Close the quantity templates
    call h5gClose_F ( gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close templates group for vector ' // trim(name) )

    ! Now setup the result
    vector = CreateVector ( GetStringIndexFromString ( trim(name) ), &
      & vt, quantities )

    ! Now loop over the quantities and read the results
    call h5gOpen_f ( vId, 'values', gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open values group for vector ' // trim(name) )
    ! Read in the quantity values
    do q = 1, vt%noQuantities
      call get_string ( vector%quantities(q)%template%name, word, &
        & strip=.true., noError=.true. )
      call LoadFromHDF5DS ( gId, trim(word), vector%quantities(q)%values )
    end do
    ! Close the group
    call h5gClose_f ( gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close values group for vector ' // trim(name) )

    ! Now loop over the quantities and perhaps read the masks
    call h5gOpen_f ( vId, 'masks', gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open masks group for vector ' // trim(name) )
    ! Read in the masks
    do q = 1, vt%noQuantities
      call get_string ( vector%quantities(q)%template%name, word, &
        & strip=.true., noError=.true. )
      if ( IsHDF5DSPresent ( gId, trim(word) ) ) then
        call CreateMask ( vector%quantities(q) )
        call LoadFromHDF5DS ( gId, trim(word), vector%quantities(q)%mask )
      end if
    end do

    ! Close the group
    call h5gClose_f ( gId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close values group for vector ' // trim(name) )

    ! Close the group
    call h5gClose_f ( vId, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group for vector ' // trim(name) )

  end subroutine ReadVectorFromHDF5
    
  ! ----------------------------------------------------------------------------    
    
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module VectorHDF5

! $Log$
! Revision 2.2  2003/05/19 22:06:31  pwagner
! Shortened names to Read..IndexFromHDF5Attr to comply with namelength standard
!
! Revision 2.1  2003/05/13 20:32:50  livesey
! First version
!
