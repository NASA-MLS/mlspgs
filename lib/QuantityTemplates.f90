! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
module QuantityTemplates         ! Quantities within vectors
!=============================================================================

  ! This module defines the `quantities' that make up vectors and their
  ! template information.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSCommon, only: NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning
  use Intrinsic, only: L_None
  use Output_m, only: Output
  use String_Table, only: Get_String

  implicit none
  public

  !------------------------------- RCS Ident Info ------------------------------
  character(len=130), private :: id = & 
       "$Id$"
  character(len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

  logical, parameter, private :: DEEBUG = .FALSE.           ! Usually FALSE

  ! Define some global parameters and data types.

  type QuantityTemplate_T

    ! Some administrative stuff

    integer :: Name            ! Sub-rosa index of quantity name
    integer :: id              ! Id code for quantity (for checking stuff)

    ! This integer is of an enumerated type describing what kind of
    ! quantity this is -- one of the l_lits of type t_quantityType
    ! in Init_Tables_Module, e.g. l_Temperature.

    integer :: quantityType

    ! The dimensions of this quantity

    integer :: noInstances     ! Number of horizontal instances in this quantity
    integer :: noSurfs         ! Number of surfaces per instance
    integer :: noChans         ! Number of channels

    ! Flags describing the quantity

    logical :: coherent        ! Do instances have same vertical coordinates?
    logical :: stacked         ! Are instances true vertical profiles?
    logical :: regular         ! Are all channels/heights represented

    ! This next one allows software using the vector quantities to be somewhat
    ! lazy and, for example, avoid interpolation.  Minor frame quantities are
    ! incoherent and unstacked, but may be regular or irregular.  However, not
    ! all incoherent unstacked quantities are minor frame quantities.

    logical :: minorFrame      ! Is this a minor frame quantity.

   ! At least in the beginning
   ! there will only be a few major frame quantities
   ! (These are vector quantities with no vert. coord. that share
   !  their other geoloc., e.g. lat and lon, with minor frame quants.)
    logical :: majorFrame      ! Is this a major frame quantity.

    ! This one indicates whether log or linear interpolation should be used
    logical :: logBasis                 ! If set use log
    real(r8) :: minValue                ! Minimum value to consider if using log

    ! This information describes how much of the data is in the overlap
    ! regions if any.

    integer :: noInstancesLowerOverlap
    integer :: noInstancesUpperOverlap

    ! Vertical coordinate

    integer :: verticalCoordinate ! The vertical coordinate used.  These
                                  ! are l_lits of the type t_VGridCoord
                                  ! defined in Init_Tables_Module.

    ! Misc. information

    real(r8) :: badValue      ! Value used to flag bad/missing data
    integer :: unit           ! Unit quantity is in when scaled as below,
                              ! an l_lit of the type t_units in Units.f90.
    real(r8) :: scaleFactor   ! Scale factor used when printing etc.

    ! For regular quantities the number of elements of each instance
    ! is simply noSurfs*noChans.  For irregular ones it is less, but it is
    ! constant from instance to instance; this is that number.

    integer :: instanceLen

    ! Give the vertical coordinates

    real(r8), dimension(:,:), pointer :: surfs => NULL()

    ! This is dimensioned (noSurfs,1) for coherent quantities and
    ! (noSurfs, noInstances) for incoherent ones.

    ! Horizontal coordinates

    real(r8), dimension(:,:), pointer :: phi => NULL()

    ! This is dimensioned (1, noInstances) for stacked quantities and
    ! (noSurfs, noInstances) for unstacked ones.
    
    ! These other coordinates are dimensioned in the same manner:

    real(r8), dimension(:,:), pointer :: geodLat => NULL()
    real(r8), dimension(:,:), pointer :: lon => NULL()
    real(r8), dimension(:,:), pointer :: time => NULL()
    real(r8), dimension(:,:), pointer :: solarTime => NULL()
    real(r8), dimension(:,:), pointer :: solarZenith => NULL()
    real(r8), dimension(:,:), pointer :: losAngle => NULL()

    ! These optional integer arrays are used for minor frame quantities,
    ! to index the major frames.

    integer, dimension(:), pointer :: mafIndex => NULL() ! (noInstances) index in l1b file
    integer, dimension(:), pointer :: mafCounter => NULL() ! (noInstances) from l1b

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! For quantities containing `channels' the following information may or
    ! may not be useful.

    ! Some quantities are on abritrary freqency grids; these quantities refer
    ! to those.

    integer :: frequencyCoordinate ! An enumerated type, e.g. FG_USBFreq
    real(r8), dimension(:), pointer :: frequencies => NULL() ! List of frequencies
                                                   ! (noChans)

    real(r8) :: lo     ! Local oscillator (optional)

    integer :: signal                   ! Index into signals database
    integer :: sideband                 ! Associated sideband -1, 0, +1

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Some families of quantities require special additional information.
    ! This is given here if needed.

    integer :: instrumentModule ! Index in the Modules database in MLSSignals_m
    integer :: radiometer       ! For ptan etc., index into radiometers database
    integer :: molecule ! What molecule does this refer to? (One of the l_...
                        ! lits of type t_molecule in Molecules.)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! For irregular quantities, instead of using the we have these arrays to
    ! help us navigate around the quantity.

    integer, dimension(:,:), pointer :: surfIndex => NULL()
    integer, dimension(:,:), pointer :: chanIndex => NULL()
    ! These are actually dimensioned (instanceLen, noInstances)
  end type QuantityTemplate_T

  ! Incrementing counter used to set the id field of a quantity template:

  integer, save, public :: quantityTemplateCounter = 0

  ! Local procedures
  interface CheckIntegrity
    module procedure CheckIntegrity_QuantityTemplate
  end interface

  private :: CheckIntegrity_QuantityTemplate
  
contains ! =====     Public Procedures     =============================

  ! Subroutines to deal with these quantitites

  ! ------------------------------  AddQuantityTemplateToDatabase  -----
  integer function AddQuantityTemplateToDatabase ( database, item )

  ! Add a quantity template to a database, or create the database if it
  ! doesn't yet exist

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: database
    type (QuantityTemplate_T), intent(in) :: item

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddQuantityTemplateToDatabase = newSize
  end function AddQuantityTemplateToDatabase

  ! ---------------------------- CheckIntegrity_QuantityTemplate -------
  logical function CheckIntegrity_QuantityTemplate ( qty, noError )
    type (QuantityTemplate_T), intent(in) :: QTY
    logical, intent(in), optional :: NOERROR
    
    ! Local variables
    integer :: NOINSTANCESOR1           ! Test value
    integer :: NOSURFSOR1               ! Test value

    integer :: MESSAGETYPE
    character ( len=132 ) :: NAME
    
    ! Executable code
    messageType = MLSMSG_Error
    if ( present ( noError ) ) then
      if ( noError ) messageType = MLSMSG_Warning
    end if

    ! Now check the integrity of the template
    if ( qty%coherent ) then
      noInstancesOr1 = 1
    else
      noInstancesOr1 = qty%noInstances
    end if

    if ( qty%stacked ) then
      noSurfsOr1 = 1
    else
      noSurfsOr1 = qty%noSurfs
    end if

    if ( qty%name /= 0 ) then
      call get_string ( qty%name, name, strip=.true. )
    else
      name = '<no name>'
    end if

    CheckIntegrity_QuantityTemplate = .true.

    ! Check the instances / overlap stuff
    if ( qty%noInstances < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noInstances for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( ( qty%noInstancesLowerOverlap < 0 ) .or. &
      &  ( qty%noInstancesLowerOverlap > qty%noInstances ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Inappropriate noInstancesLowerOverlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( ( qty%noInstancesUpperOverlap < 0 ) .or. &
      &  ( qty%noInstancesUpperOverlap > qty%noInstances ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Inappropriate noInstancesUpperOverlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    if ( qty%noInstancesLowerOverlap + qty%noInstancesUpperOverlap > &
      & qty%noInstances ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Too much overlap for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check the surfaces stuff
    if ( qty%noSurfs < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noSurfs for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check the channels stuff
    if ( qty%noChans < 0 ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad noChans for quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check the instanceLen
    if ( qty%regular .and. &
      & qty%InstanceLen /= qty%noSurfs * qty%noChans ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The quantity template '//trim(name)//' does not have the right instanceLen' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Check the arrays are associated.  Note these have to be errors, as later
    ! tests will fail otherwise.
    if ( .not. associated ( qty%surfs ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have surfs associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( .not. associated ( qty%phi ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have phi associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( .not. associated ( qty%geodLat ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have geodLat associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( .not. associated ( qty%lon ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have lon associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( .not. associated ( qty%time ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have time associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( .not. associated ( qty%solarTime ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have solarTime associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( .not. associated ( qty%solarZenith ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have solarZenith associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( .not. associated ( qty%losAngle ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity template '//trim(name)// ' does not have losAngle associated' )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check the array lower bounds
    if ( any ( lbound ( qty%surfs ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for surfs array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( any ( lbound ( qty%phi ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for phi array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( any ( lbound ( qty%geodLat ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for geodLat array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( lbound ( qty%lon ) /= (/1,1/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for lon array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( lbound ( qty%time ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for time array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( lbound ( qty%solarTime ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for solarTime array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( lbound ( qty%solarZenith ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for solarZenith array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( lbound ( qty%losAngle ) /= 1 ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad lbound for losAngle array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check the array upper bounds
    if ( any ( ubound ( qty%surfs ) /= (/qty%noSurfs, noInstancesOr1/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for surfs array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( any ( ubound ( qty%phi ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for phi array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    if ( any ( ubound ( qty%geodLat ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for geodLat array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( ubound ( qty%lon ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for lon array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( ubound ( qty%time ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for time array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( ubound ( qty%solarTime ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for solarTime array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( ubound ( qty%solarZenith ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for solarZenith array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if
    if ( any ( ubound ( qty%losAngle ) /= (/noSurfsOr1, qty%noInstances/) ) ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'Bad ubound for losAngle array in quantity template '//trim(name) )
      CheckIntegrity_QuantityTemplate = .false.      
    end if

    ! Check stuff specific for major/minor frame stuff
    if ( qty%minorFrame .or. qty%majorFrame ) then
      if ( .not. associated ( qty%mafIndex ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'The quantity template '//trim(name)// ' does not have mafIndex associated' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if
      if ( .not. associated ( qty%mafCounter ) ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'The quantity template '//trim(name)// ' does not have mafCounter associated' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if

      if ( lbound ( qty%mafIndex, 1 ) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Bad lbound for mafIndex' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if
      if ( lbound ( qty%mafCounter, 1 ) /= 1 ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Bad lbound for mafCounter' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if

      if ( ubound ( qty%mafIndex, 1 ) /= qty%noInstances ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Bad ubound for mafIndex' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if
      if ( ubound ( qty%mafIndex, 1 ) /= qty%noInstances ) then
        call MLSMessage ( messageType, ModuleName, &
          & 'Bad ubound for mafIndex' )
        CheckIntegrity_QuantityTemplate = .false.      
      end if
    end if

    ! Check irregular stuff
    if ( .not. qty%regular ) then
      call MLSMessage ( messageType, ModuleName, &
        & 'The quantity '//trim(name)//' appears to be irregular' )
      CheckIntegrity_QuantityTemplate = .false.
    end if

    ! Could check channels stuff here, but not sure what to do.

  end function CheckIntegrity_QuantityTemplate

  ! ----------------------------  DestroyQuantityTemplateContents  -----
  subroutine DestroyQuantityTemplateContents ( qty )

  ! Destroy a quantity template

    ! Dummy argument
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Local variables
    character (LEN=32) :: quantityNameStr

    ! Executable code

    if ( DEEBUG ) then
      call output('Destroying Quantity template contents', advance='yes')
      call output('minor frame? ', advance='no')
      call output(qty%minorFrame, advance='no')
      call output('   major frame? ', advance='no')
      call output(qty%majorFrame, advance='no')
      call output('   template  name ', advance='no')
      if ( qty%name < 1 ) then
        call output('   (unnamed) ', advance='yes')
      else
        call get_string(qty%name, quantityNameStr, strip=.true., noerror=.true.)
        call output(trim(quantityNameStr), advance='yes')
      endif
    endif
    if ( DEEBUG ) call output('Deallocating qty%surfs', advance='no')
    call deallocate_test ( qty%surfs, "qty%surfs", ModuleName )
    if ( DEEBUG ) call output('  qty%phi', advance='no')
    call deallocate_test ( qty%phi, "qty%phi", ModuleName )
    if ( DEEBUG ) call output('  qty%geodLat', advance='no')
    call deallocate_test ( qty%geodLat, "qty%geodLat", ModuleName )
    if ( DEEBUG ) call output('  qty%lon', advance='no')
    call deallocate_test ( qty%lon, "qty%lon", ModuleName )
    if ( DEEBUG ) call output('  qty%time', advance='no')
    call deallocate_test ( qty%time, "qty%time", ModuleName )
    if ( DEEBUG ) call output('  qty%solartime', advance='no')
    call deallocate_test ( qty%solarTime, "qty%solarTime", ModuleName )
    if ( DEEBUG ) call output('  qty%solarzenith', advance='no')
    call deallocate_test ( qty%solarZenith, "qty%solarZenith", ModuleName )
    if ( DEEBUG ) call output('  qty%losAngle', advance='no')
    call deallocate_test ( qty%losAngle, "qty%losAngle", ModuleName )
    if ( DEEBUG ) call output('  qty%frequencies', advance='yes')
    call deallocate_test ( qty%frequencies, "qty%frequencies", ModuleName )

!    if (qty%minorFrame .or. qty%majorFrame) then
    if ( qty%minorFrame ) then
      if ( DEEBUG ) call output('Deallocating qty%MAFIndex', advance='no')
      call deallocate_test ( qty%MAFIndex, "qty%MAFIndex", ModuleName )
      if ( DEEBUG ) call output('  qty%MAFCounter', advance='yes')
      call deallocate_test ( qty%MAFCounter, "qty%MAFCounter", ModuleName )
    end if
    
    if (.NOT. qty%regular) then
      if ( DEEBUG ) call output('Deallocating qty%surfIndex', advance='no')
      call deallocate_test ( qty%surfIndex, "qty%surfIndex", ModuleName )
      if ( DEEBUG ) call output('  qty%chanIndex', advance='yes')
      call deallocate_test ( qty%chanIndex, "qty%chanIndex", ModuleName )
    end if

  end subroutine DestroyQuantityTemplateContents

  ! ----------------------------  DestroyQuantityTemplateDatabase  -----
  subroutine DestroyQuantityTemplateDatabase ( database, &
    & ignoreMinorFrame, ignoreMajorFrame )

  ! Destroy a quantity template database

    ! Dummy argument
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE
    logical, intent(in), optional :: ignoreMinorFrame
    logical, intent(in), optional :: ignoreMajorFrame

    ! Local variables
    integer :: qtyIndex, status
    logical :: myIgnoreMinorFrame
    logical :: myIgnoreMajorFrame
    
    myIgnoreMinorFrame = .false.
    if ( present ( ignoreMinorFrame ) ) myIgnoreMinorFrame = ignoreMinorFrame
    myIgnoreMajorFrame = .false.
    if ( present ( ignoreMajorFrame ) ) myIgnoreMajorFrame = ignoreMajorFrame

    if ( associated(database) ) then
      do qtyIndex = 1, SIZE(database)
        if ( &
          & .not. ( &
            & (database(qtyIndex)%minorFrame .and. myIgnoreMinorFrame) &
            & .or. &
            & (database(qtyIndex)%majorFrame .and. myIgnoreMajorFrame) &
            & ) &
          & ) call DestroyQuantityTemplateContents ( database(qtyIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if
  end subroutine DestroyQuantityTemplateDatabase

  ! ---------------------------------- InflateQuantityTemplateDatabase --
  integer function InflateQuantityTemplateDatabase ( database, extra )
    ! Make a quantity template database bigger by extra
    ! Return index of first new element

    ! Dummy arguments
    type (QuantityTemplate_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type (QuantityTemplate_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateQuantityTemplateDatabase = firstNewItem
  end function InflateQuantityTemplateDatabase

  ! -----------------------------------  SetupNewQuantityTemplate  -----
  subroutine SetupNewQuantityTemplate ( qty, source, noInstances, noSurfs, &
    & noChans, coherent, stacked, regular, instanceLen, minorFrame, majorFrame )

  ! Set up a new quantity template according to the user input.  This may
  ! be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

    ! Dummy arguments
    type (QuantityTemplate_T), intent(inout) :: qty ! Result

    type (QuantityTemplate_T), optional, intent(in) :: source ! Template
    integer, intent(in), optional :: noInstances
    integer, intent(in), optional :: noSurfs
    integer, intent(in), optional :: noChans
    logical, intent(in), optional :: coherent
    logical, intent(in), optional :: stacked
    logical, intent(in), optional :: regular
    integer, intent(in), optional :: instanceLen
    logical, intent(in), optional :: minorFrame
    logical, intent(in), optional :: majorFrame

    ! Local variables
    integer :: noSurfsToAllocate        ! For allocations
    integer :: noInstancesToAllocate    ! For allocations

    ! Executable code

    call destroyQuantityTemplateContents ( qty ) ! Avoid memory leaks

    ! First, if we have a template setup according to that
    if (present(source)) then
      qty%noInstances = source%noInstances
      qty%noSurfs = source%noSurfs
      qty%noChans = source%noChans
      qty%coherent = source%coherent
      qty%stacked = source%stacked
      qty%regular = source%regular
      qty%minorFrame = source%minorFrame
      qty%majorFrame = source%majorFrame
      qty%logBasis = source%logBasis
      qty%instanceLen =  source%instanceLen
      qty%verticalCoordinate = source%verticalCoordinate
      qty%frequencyCoordinate = source%frequencyCoordinate
    else ! We have no template, setup a very bare quantity
      qty%noInstances = 1
      qty%noSurfs = 1
      qty%noChans = 1
      qty%coherent = .TRUE.
      qty%stacked = .TRUE.
      qty%regular = .TRUE.
      qty%logBasis = .FALSE.
      qty%minorFrame = .FALSE.
      qty%majorFrame = .FALSE.
      qty%instanceLen = 1
      qty%verticalCoordinate=l_none
      qty%frequencyCoordinate=l_none
    end if

    ! Now, see if the user asked for modifications to this
    if ( present(noInstances) ) qty%noInstances = noInstances
    if ( present(noSurfs) ) qty%noSurfs = noSurfs
    if ( present(noChans) ) qty%noChans = noChans
    if ( present(regular) ) qty%regular = regular
    if ( present(minorFrame) ) qty%minorFrame = minorFrame
    if ( present(majorFrame) ) qty%majorFrame = majorFrame
    if ( qty%minorFrame ) then
      if ( present(coherent) ) then
        if ( coherent ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Minor frame quantities must be incoherent" )
      end if
      qty%coherent = .FALSE.
      if ( present(stacked) ) then
        if ( stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Minor frame quantities must be unstacked" )
      end if
      qty%stacked = .FALSE.
    else
      if ( present(coherent) ) qty%coherent = coherent
      if ( present(stacked) ) qty%stacked = stacked
    end if

    ! Now think about instanceLen
    if ( (.NOT. qty%regular) .AND. (present(instanceLen)) ) then
      qty%instanceLen = instanceLen
    else
      qty%instanceLen = qty%noSurfs*qty%noChans
    end if

    ! Now we allocate all the arrays we're going to need

    if ( qty%coherent ) then 
      noInstancesToAllocate = 1
    else
      noInstancesToAllocate = qty%noInstances
    end if

    if ( qty%stacked ) then
      noSurfsToAllocate = 1
    else
      noSurfsToAllocate = qty%noSurfs
    end if

    ! First the vertical coordinates

    call allocate_test ( qty%surfs, qty%noSurfs, noInstancesToAllocate, &
      & "qty%surfs", ModuleName )

    ! Now the horizontal coordinates

    call allocate_test ( qty%phi, noSurfsToAllocate, qty%noInstances, &
      & "qty%phi", ModuleName )

    call allocate_test ( qty%geodLat, noSurfsToAllocate, qty%noInstances, &
      & "qty%geodLat", ModuleName )

    call allocate_test ( qty%lon, noSurfsToAllocate, qty%noInstances, &
      & "qty%lon", ModuleName )

    call allocate_test ( qty%time, noSurfsToAllocate, qty%noInstances, &
      & "qty%time", ModuleName )

    call allocate_test ( qty%solarTime, noSurfsToAllocate, qty%noInstances, &
      & "qty%solarTime", ModuleName )

    call allocate_test ( qty%solarZenith, noSurfsToAllocate, qty%noInstances, &
      & "qty%solarZenith", ModuleName )

    call allocate_test ( qty%losAngle, noSurfsToAllocate, qty%noInstances, &
      & "qty%losAngle", ModuleName )

    ! Now some other stuff to allocate

    if ( qty%minorFrame ) then
      call allocate_test ( qty%MAFIndex, qty%noInstances, &
        & "qty%MAFIndex", ModuleName )
      call allocate_test ( qty%MAFCounter, qty%noInstances, &
        & "qty%MAFCounter", ModuleName )
    else
      nullify ( qty%MAFIndex, qty%MAFCounter )
    end if

    if (.NOT. qty%regular) then
      call allocate_test ( qty%surfIndex, qty%instanceLen, qty%noInstances, &
        & "qty%surfIndex", ModuleName )
      call allocate_test ( qty%chanIndex, qty%instanceLen, qty%noInstances, &
        & "qty%chanIndex", ModuleName )
    else
      nullify ( qty%surfIndex, qty%chanIndex )
    end if

    ! Increment the id counter and set the id field
    quantityTemplateCounter = quantityTemplateCounter + 1
    qty%id = quantityTemplateCounter
  end subroutine SetupNewQuantityTemplate

!=============================================================================
end module QuantityTemplates
!=============================================================================

!
! $Log$
! Revision 2.25  2002/08/28 20:42:11  livesey
! Added InflateQuantityTemplateDatabase
!
! Revision 2.24  2002/07/22 03:26:05  livesey
! Added checkIntegrity
!
! Revision 2.23  2002/07/01 23:51:07  vsnyder
! Plug a memory leak
!
! Revision 2.22  2001/10/12 23:09:25  pwagner
! More debugging statements
!
! Revision 2.21  2001/10/03 17:42:27  pwagner
! reset DEEBUG to FALSE
!
! Revision 2.20  2001/10/02 23:12:50  pwagner
! More chi^2 fixes
!
! Revision 2.19  2001/09/17 21:59:26  livesey
! Removed allocate of frequencies, it's deferred to later in the code
!
! Revision 2.18  2001/09/13 19:59:43  pwagner
! Added majorframe as possible quantity type
!
! Revision 2.17  2001/07/31 23:39:12  dwu
! allocate and deallocate qty%frequencies
!
! Revision 2.16  2001/07/11 21:41:16  livesey
! Made quantityTemplateCounter public
!
! Revision 2.15  2001/07/02 17:25:30  livesey
! Some changes to comments, following walk through
!
! Revision 2.14  2001/05/23 20:38:35  livesey
! Updated a comment
!
! Revision 2.13  2001/04/23 23:52:16  livesey
! Sorry, should have put comment in one below.  Now has optional ignoreMinorFrame
! argument to DestroyQuantityTemplateDatabase
!
! Revision 2.12  2001/04/23 23:50:41  livesey
! *** empty log message ***
!
! Revision 2.11  2001/04/12 21:43:06  livesey
! Added sideband field
!
! Revision 2.10  2001/04/10 22:37:49  vsnyder
! Fix a type
!
! Revision 2.9  2001/03/24 00:31:12  pwagner
! USEs output in case we replace MLSMessage with output in additem..
!
! Revision 2.8  2001/03/17 02:23:18  livesey
! Added log basis field
!
! Revision 2.7  2001/03/15 20:20:59  vsnyder
! Correct the description of 'InstrumentModule'
!
! Revision 2.6  2001/03/02 01:34:03  livesey
! New signals stuff
!
! Revision 2.5  2001/02/23 17:47:01  livesey
! Nullified pointers.
!
! Revision 2.4  2001/02/14 00:12:34  livesey
! Removed firstIndexChannel
!
! Revision 2.3  2001/02/09 00:38:56  livesey
! Various changes
!
! Revision 2.2  2000/12/04 23:43:59  vsnyder
! Move more of addItemToDatabase into the include
!
! Revision 2.1  2000/10/13 00:00:37  vsnyder
! Moved from mlspgs/l2 to mlspgs/lib
!
! Revision 2.0  2000/09/05 18:57:04  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
