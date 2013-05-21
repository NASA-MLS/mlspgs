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
module FillUtils_1                     ! Procedures used by Fill
  !=============================================================================

  use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use CHUNKS_M, only: MLSCHUNK_T
  use CONSTANTS, only: DEG2RAD, LN10, RAD2DEG
  use DUMP_0, only: DUMP
  use EXPR_M, only: EXPR, EXPR_CHECK, GETINDEXFLAGSFROMLIST
  use GRIDDEDDATA, only: GRIDDEDDATA_T, DUMP, WRAPGRIDDEDDATA
  use HDF5, only: HSIZE_T, H5DOPEN_F, H5DCLOSE_F
  use INIT_TABLES_MODULE, only: F_MEASUREMENTS, F_TOTALPOWERVECTOR, &
    & F_WEIGHTSVECTOR, &
    & L_ADDNOISE, L_BASELINE, L_BINMAX, L_BINMEAN, L_BINMIN, L_BINTOTAL, &
    & L_BOUNDARYPRESSURE, L_BOXCAR, L_CHISQBINNED, L_CHISQCHAN, L_CHISQRATIO, &
    & L_CHISQMMAF, L_CHISQMMIF, &
    & L_CLOUDINDUCEDRADIANCE, L_CLOUDMINMAX, &
    & L_COLUMNABUNDANCE, &
    & L_DNWT_FLAG, L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, L_DNWT_CHISQRATIO, &
    & L_DOBSONUNITS, L_DU, &
    & L_ECRTOFOV, &
    & L_FIELDAZIMUTH, L_FIELDELEVATION, L_FIELDSTRENGTH, &
    & L_GPH, &
    & L_HEIGHT, L_ISOTOPERATIO, &
    & L_L1BMAFBASELINE, L_L1BMIF_TAI, &
    & L_LIMBSIDEBANDFRACTION, L_LOSVEL, &
    & L_LSLOCAL, L_LSGLOBAL, L_LSWEIGHTED, &
    & L_MAGNETICFIELD, &
    & L_MAX, L_MEAN, L_MIN, L_MOLCM2, &
    & L_NOISEBANDWIDTH, L_NONE, L_NORADSPERMIF, L_NORADSBINNED, &
    & L_ORBITINCLINATION, &
    & L_PRESSURE, L_PTAN,  L_QUALITY, &
    & L_RADIANCE, L_REFGPH, &
    & L_REFLTEMP, &
    & L_SCECI, L_SCGEOCALT, L_SCVEL, L_SCVELECI, L_SCVELECR, &
    & L_SINGLECHANNELRADIANCE, &
    & L_STATUS, L_SYSTEMTEMPERATURE, &
    & L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, &
    & L_TNGTGEOCALT, L_TOTALPOWERWEIGHT, L_VMR, &
    & L_XYZ, L_ZETA
  use INTRINSIC, only: FIELD_INDICES, LIT_INDICES, &
    & PHYQ_DIMENSIONLESS, PHYQ_INDICES, PHYQ_INVALID, PHYQ_TEMPERATURE, &
    & PHYQ_LENGTH, PHYQ_PRESSURE, PHYQ_ZETA, PHYQ_ANGLE
  use L1BDATA, only: DEALLOCATEL1BDATA, DUMP, GETL1BFILE, L1BDATA_T, &
    & PRECISIONSUFFIX, READL1BDATA, ASSEMBLEL1BQTYNAME
  use L2GPDATA, only: L2GPDATA_T, READl2GPDATA, DESTROYl2GPCONTENTS
  use L2AUXDATA, only: L2AUXDATA_T, MAXSDNAMESBUFSIZE, &
    & READL2AUXDATA, DESTROYL2AUXCONTENTS
  use L3ASCII, only: L3ASCII_INTERP_FIELD
  use MANIPULATEVECTORQUANTITIES, only: DOFGRIDSMATCH, DOHGRIDSMATCH, &
    & DOVGRIDSMATCH, DOQTYSDESCRIBESAMETHING
  use MATRIXMODULE_0, only: MATRIXELEMENT_T, M_FULL, &
    & CREATEBLOCK, SPARSIFY, MATRIXINVERSION
  use MATRIXMODULE_1, only: DUMP, FINDBLOCK, MATRIX_SPD_T, UPDATEDIAGONAL
  ! NOTE: IF YOU EVER WANT TO INCLUDE DEFINED ASSIGNMENT FOR MATRICES, PLEASE
  ! CAREFULLY CHECK OUT THE CODE AROUND THE CALL TO SNOOP.
  use MLSCOMMON, only: MLSFILE_T, DEFAULTUNDEFINEDVALUE
  use MLSFILES, only: HDFVERSION_5, DUMP, GETMLSFILEBYTYPE
  use MLSFILLVALUES, only: ISFILLVALUE, ISMONOTONIC, MONOTONIZE, REMOVEFILLVALUES
  use MLSKINDS, only: R4, R8, RM, RP, RV
  use MLSL2OPTIONS, only: L2CFNODE, MLSMESSAGE
  use MLSMESSAGEMODULE, only: MLSMSG_ERROR, MLSMSG_WARNING, &
    & MLSMESSAGECALLS
  use MLSNUMERICS, only: COEFFICIENTS_R8, INTERPOLATEARRAYSETUP, &
    & INTERPOLATEARRAYTEARDOWN, INTERPOLATEVALUES, HUNT
  use MLSSETS, only: FINDFIRST, FINDLAST
  use MLSSIGNALS_M, only: GETFIRSTCHANNEL, GETSIGNALNAME, GETMODULENAME, &
    & GETSIGNAL, ISMODULESPACECRAFT, SIGNAL_T, SIGNALS
  use MLSSTRINGLISTS, only: NUMSTRINGELEMENTS, &
    & STRINGELEMENT, SWITCHDETAIL
  use MLSSTRINGS, only: INDEXES, LOWERCASE, WRITEINTSTOCHARS
  use MOLECULES, only: L_H2O
  use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT, OUTPUTNAMEDVALUE
  use QUANTITYTEMPLATES, only: EPOCH, QUANTITYTEMPLATE_T
  use RHIFROMH2O, only: H2OPRECFROMRHI, RHIFROMH2O_FACTOR, RHIPRECFROMH2O
  use SCANMODELMODULE, only: GETBASISGPH, GET2DHYDROSTATICTANGENTPRESSURE, &
    & GETGPHPRECISION
  use SPECTROSCOPYCATALOG_M, only: CATALOG
  use STRING_TABLE, only: DISPLAY_STRING, GET_STRING
  use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, SUBTREE, NSONS, &
    & SOURCE_REF, SUBTREE
  use VECTORSMODULE, only: &
    & CLEARUNDERMASK, CLONEVECTORQUANTITY, COPYVECTOR, CREATEMASK, &
    & DESTROYVECTORINFO, DESTROYVECTORQUANTITYMASK, &
    & DESTROYVECTORQUANTITYVALUE, DUMP, &
    & GETVECTORQTYBYTEMPLATEINDEX, GETVECTORQUANTITYBYTYPE, &
    & ISVECTORQTYMASKED, MASKVECTORQTY, &
    & VALIDATEVECTORQUANTITY, VECTOR_T, &
    & VECTORVALUE_T, M_CLOUD, M_FILL, M_IGNORE, M_LINALG
  use VGRIDSDATABASE, only: VGRID_T, GETUNITFORVERTICALCOORDINATE

  implicit none
  private
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  logical, parameter :: COUNTEMPTY = .true.
  logical, parameter :: DONTPAD = .false.
  logical, parameter :: WARNWHENPTANNONMONOTONIC = .false.

  ! -----     Declarations for Fill and internal subroutines     -------

  logical, parameter :: ADDSLASH = .false.
  logical, parameter :: DEEBUG = .FALSE.                 ! Usually FALSE
  logical, parameter :: UNIFORMCHISQRATIO = .FALSE.
  integer, public :: FILLERROR

  ! -999.99 ! Same as %template%badvalue
  real, parameter ::    UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE

  ! Error codes for "announce_error"
  integer, parameter, public :: No_Error_code = 0
  integer, parameter, public :: CantFromL2AUX = No_Error_code + 1
  integer, parameter, public :: CantFromL1B = cantFromL2AUX + 1

  ! Error codes for "Matrix" specification
  integer, parameter, public :: MissingField = cantFromL1B + 1

  ! More Error codes relating to Vector
  integer, parameter, public :: NumChansisZero = missingField + 1
  integer, parameter, public :: NoSourceGridGiven= numChansisZero + 1
  integer, parameter, public :: NoSourceL2GPGiven= noSourceGridGiven + 1
  integer, parameter, public :: NoSourceL2AUXGiven= noSourceL2GPGiven + 1
  integer, parameter, public :: NoExplicitValuesGiven= noSourceL2AUXGiven + 1
  integer, parameter, public :: InvalidExplicitFill = noExplicitValuesGiven + 1
  integer, parameter, public :: BadIsotopeFill = invalidExplicitFill + 1
  integer, parameter, public :: BadlosGridFill = badIsotopeFill + 1
  integer, parameter, public :: CantInterpolate3d = badlosGridFill + 1
  integer, parameter, public :: WrongUnits = CantInterpolate3d + 1

  ! Error codes resulting from Covariance
  integer, parameter, public :: BothFractionAndLength = WrongUnits + 1
  integer, parameter, public :: Inconsistent = BothFractionAndLength + 1
  integer, parameter, public :: NotImplemented = Inconsistent + 1
  integer, parameter, public :: NotPlain = NotImplemented  + 1
  integer, parameter, public :: NotSPD = NotPlain + 1

  ! Miscellaneous
  integer, parameter, public :: BadGeocAltitudeQuantity = NotSPD + 1
  integer, parameter, public :: BadTemperatureQuantity = badGeocAltitudeQuantity + 1
  integer, parameter, public :: BadREFGPHQuantity = badTemperatureQuantity + 1
  integer, parameter, public :: Miscellaneous_err = badREFGPHQuantity + 1
  integer, parameter, public :: ErrorReadingL1B = miscellaneous_err + 1
  integer, parameter, public :: NeedTempREFGPH = errorReadingL1B + 1
  integer, parameter, public :: NeedH2O = needTempRefGPH + 1
  integer, parameter, public :: NeedOrbitInclination = needH2O + 1
  integer, parameter, public :: NeedGeocAltitude = needOrbitInclination + 1
  integer, parameter, public :: NoCodeFor = needGeocAltitude + 1
  integer, parameter, public :: NonConformingHydrostatic = noCodeFor + 1
  integer, parameter, public :: NoSpecialFill = nonConformingHydrostatic + 1
  integer, parameter, public :: NoTotalPower = noSpecialFill + 1
  integer, parameter, public :: BadlosVelFill = noTotalPower + 1
  integer, parameter, public :: NotZetaForGrid = BadLosVelFill + 1
  integer, parameter, public :: BadEstNoiseFill = NotZetaForGrid + 1
  integer, parameter, public :: BadRefractFill = BadEstNoiseFill + 1
  integer, parameter, public :: MissingDataInGrid = BadRefractFill + 1
  integer, parameter, public :: EmptyGridForFill = MissingDataInGrid + 1

  logical :: UNITSERROR               ! From expr

  public :: ADDGAUSSIANNOISE, APPLYBASELINE, AUTOFILLVECTOR, &
      & COMPUTETOTALPOWER, DEALLOCATESTUFF, &
      & EXTRACTSINGLECHANNEL, FILLCOVARIANCE, FROMANOTHER, FROMGRID, &
      & FROML2GP, FROMPROFILE, LOSVELOCITY, &
      & CHISQCHAN, CHISQMMAF, CHISQMMIF, CHISQRATIO, &
      & COLABUNDANCE, DERIVATIVEOFSOURCE, FOLDEDRADIANCE, PHITANWITHREFRACTION, &
      & IWCFROMEXTINCTION, RHIFROMORTOH2O, NORADSPERMIF, &
      & RHIPRECISIONFROMORTOH2O, WITHESTNOISE, &
      & HYDROSTATICALLY, FROMSPLITSIDEBAND, GPHPRECISION, &
      & FROMISOTOPE, FROMASCIIFILE, ROTATEMAGNETICFIELD, &
      & EXPLICIT, FROML1B, &
      & FROML2AUX, USINGMAGNETICMODEL, &
      & FROMINTERPOLATEDQTY, FROMLOSGRID, &
      & BYMANIPULATION, MANIPULATEVECTORS, WITHREFLECTORTEMPERATURE, &
      & WITHREICHLERWMOTP, &
      & WITHWMOTROPOPAUSE, WITHBINRESULTS, WITHBOXCARFUNCTION, &
      & STATUSQUANTITY, QUALITYFROMCHISQ, CONVERGENCEFROMCHISQ, &
      & USINGLEASTSQUARES, OFFSETRADIANCEQUANTITY, RESETUNUSEDRADIANCES, &
      & SCALEOVERLAPS, SPREADCHANNELFILL, &
      & TRANSFERVECTORS, TRANSFERVECTORSBYMETHOD, &
      & UNCOMPRESSRADIANCE, &
      & ANNOUNCE_ERROR, QTYFROMFILE, VECTORFROMFILE

  interface FromProfile
    module procedure FromProfile_node, FromProfile_values
  end interface

  logical, parameter :: WARNIFVERTCOORDNOTZETA = .false.

contains ! =====     Public Procedures     =============================

    ! ------------------------------------------- addGaussianNoise ---
    subroutine addGaussianNoise ( key, quantity, sourceQuantity, &
              & noiseQty, multiplier, spread, ignoreTemplate )
      use MLSRandomNumber, only: DRANG
      ! A special fill: quantity = sourceQuantity + g() noiseQty
      ! where g() is a random number generator with mean 0 and std. dev. 1
      ! Generalized into ( a sourceQuantity + b g() noiseQty )
      ! where a and b are multipliers)
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) ::      quantity
      type (VectorValue_T), intent(in) ::         sourceQuantity
      type (VectorValue_T), intent(in) ::         noiseQty
      real, dimension(:), intent(in) ::           multiplier
      logical, intent(in)              ::         ignoreTemplate
      logical, intent(in)              ::         spread

      ! Local variables
      integer                          ::    ROW, COLUMN
      real                             ::    a, b

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.addGaussianNoise', key )
      ! First check that things are OK.
      if ( .not. ignoreTemplate .and. .not. FillableChiSq ( quantity, &
        & sourceQuantity, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities adding noise'  )
        go to 9
      end if

     ! Either multiplier = [a, b] or multiplier = b are possible
      if ( &
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 1.
        b = multiplier(1)
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      if ( spread .and. &
        & size(quantity%values, 1) /= size(sourceQuantity%values, 1) ) then
        do column=1, size(quantity%values(1, :))
          quantity%values(:, column) = &
            & sourceQuantity%values(1, column) * a &
            & + &
            & drang() * noiseQty%values(1, column) * b
        end do
      elseif ( spread .and. &
        & size(quantity%values, 2) /= size(sourceQuantity%values, 2) ) then
        do row=1, size(quantity%values(:, 1))
          quantity%values(row, :) = &
            & sourceQuantity%values(row, :) * a &
            & + &
            & drang() * noiseQty%values(row, :) * b
        end do
      else
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            quantity%values(row, column) = &
              & sourceQuantity%values(row, column) * a &
              & + &
              & drang() * noiseQty%values(row, column) * b
          end do
        end do
      endif

    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.addGaussianNoise' )
    end subroutine addGaussianNoise

    ! ---------------------------------------------  ANNOUNCE_ERROR  -----
    subroutine ANNOUNCE_ERROR ( WHERE, CODE, EXTRAMESSAGE, EXTRAINFO, QUITNOW )

      use MORETREE, only: GET_FIELD_ID, STARTERRORMESSAGE

      integer, intent(in) :: WHERE   ! Tree node where error was noticed
      integer, intent(in) :: CODE    ! Code for error message
      character (len=*), intent(in), optional :: EXTRAMESSAGE
      integer, intent(in), dimension(:), optional :: EXTRAINFO
      logical, intent(in), optional :: QUITNOW

      integer :: I

      fillerror = max(fillerror,1)
      if ( present(extraMessage) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & trim(extraMessage) )
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Calling ANNOUNCE_ERROR' )
      endif
      call StartErrorMessage ( where )
      if ( code  > no_Error_Code ) call output ( 'The' );

      select case ( code )
      case ( badEstNoiseFill )
        call output ( " estimated noise fill is missing information", advance='yes' )
      case ( badGeocAltitudeQuantity )
        call output ( " geocAltitudeQuantity is not geocAltitude", advance='yes' )
      case ( badlosGridfill )
        call output ( " information for los Grid fill is incomplete/incorrect", advance='yes' )
      case ( badlosvelfill )
        call output ( " information for los velocity is incomplete/incorrect", advance='yes' )
      case ( badIsotopeFill )
        call output ( " information for isotope fill is incomplete/incorrect", advance='yes' )
      case ( badREFGPHQuantity )
        call output ( " refGPHQuantity is not refGPH", advance='yes' )
      case ( badRefractFill )
        call output ( " phiTan refract fill is missing information", advance='yes' )
      case ( badTemperatureQuantity )
        call output ( " temperatureQuantity is not temperature", advance='yes' )
      case ( bothFractionAndLength )
        call output ( " config specifies both fraction and lengthScale", advance='yes' )
      case ( cantFromL1B )
        call output ( " command could not be filled from L1B.", advance='yes' )
      case ( cantFromL2AUX )
        call output ( " command could not be filled from L2AUX.", advance='yes' )
      case ( cantInterpolate3D )
        call output ( " program cannot interpolate 3d quantities (yet).", advance='yes' )
      case ( emptyGridForFill )
        call output ( " config specifies an empty grid for the fill", advance='yes' )
      case ( errorReadingL1B )
        call output ( " L1B file could not be read.", advance='yes' )
      case ( inconsistent )
        call output ( " matrix and vector are inconsistent.", advance='yes' )
      case ( invalidExplicitFill )
        call output ( " value has inappropriate dimensionality for explicit fill.", advance='yes' )
        call output ( " Should have " )
        call output ( extraInfo )
        call output ( " elements.", advance='yes' )
      case ( missingDataInGrid )
        call output ( " grid for fill has missing/bad data points", advance='yes' )
      case ( missingField )
        call output ( " fields " )
        do i = 1, size(extraInfo)
          call display_string ( field_indices(extraInfo(i)) )
          if ( i == size(extraInfo) ) then
            call output ( " and " )
          else
            call output ( ", " )
          end if
        end do
        call output ( " are required.", advance='yes' )
      case ( noCodeFor )
        call output ( " no code for " )
        call display_string ( field_indices(extraInfo(i)), advance='yes' )
      case ( needGeocAltitude )
        call output ( " fill needs geocAltitudeQuantity.", advance='yes' )
      case ( needH2O )
        call output ( " fill needs H2OQuantity.", advance='yes' )
      case ( needOrbitInclination )
        call output ( " fill needs OrbitalInclination.", advance='yes' )
      case ( needTempREFGPH )
        call output ( " needs temperatureQuantity and refGPHquantity.", advance='yes' )
      case ( no_Error_Code ) ! Handled at the bottom
      case ( noExplicitValuesGiven )
        call output ( " explicit fill requires explicit values.", advance='yes' )
      case ( nonConformingHydrostatic )
        call output ( " quantities needed for hydrostatic fill do not conform", advance='yes' )
      case ( noSourceGridGiven )
        call output ( " gridded fill requires a sourceGrid field.", advance='yes' )
      case ( noSourceL2AUXGiven )
        call output ( " L2AUX fill requires a sourceL2AUX field.", advance='yes' )
      case ( noSourceL2GPGiven )
        call output ( " L2GP fill requires a sourceL2GP field.", advance='yes' )
      case ( noSpecialFill )
        call output ( " special fill is invalid", advance='yes' )
      case ( noTotalPower )
        call output ( " total power weights are zero", advance='yes' )
      case ( notImplemented )
        call output ( extraMessage )
        call output ( " method is not implemented yet.", advance='yes' )
      case ( notSPD )
        call output ( " matrix is not a SPD matrix.", advance='yes' )
      case ( notPlain )
        call output ( " matrix is not a plain matrix.", advance='yes' )
      case ( NotZetaForGrid )
        call output ( " Quantity not on zeta surfaces.", advance='yes' )
      case ( wrongUnits )
        call display_string ( field_indices(Get_Field_Id(where)), &
          & before=" values of the " )
        call output ( " field have the wrong units", advance="yes" )
        if ( present(extraInfo) ) then
          i = 1
          if ( size(extraInfo) > 1 ) then
            call display_String ( phyq_indices(extraInfo(1)), &
              & before="Units are " )
            call output ( ", " )
            i = 2
          end if
          call display_String ( phyq_indices(extraInfo(i:)), advance='yes', &
            & before="Units should be" )
        end if
      case default
        call outputNamedValue( ' error code', Code, advance='no' )
        call output ( " command caused an unrecognized programming error", advance='yes' )
      end select
      if ( present(ExtraMessage) )  call output(ExtraMessage, advance='yes')
      if ( code == no_Error_Code .and. present(extraInfo) ) &
        & call dump ( extraInfo, name='Extra info' )
      if ( present(QUITNOW) ) then
        if ( QUITNOW ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & trim(extraMessage) )
      end if
    end subroutine ANNOUNCE_ERROR

    ! ------------------------------------------- ApplyBaseline ----------
    subroutine ApplyBaseline ( key, quantity, baselineQuantity, &
      & quadrature, dontmask, ignoreTemplate )
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Radiance quantity to modify
      type (VectorValue_T), intent(in) :: BASELINEQUANTITY ! L1B MAF baseline to use
      logical, intent(in) :: QUADRATURE ! If set add in quadrature (for noise)
      logical, intent(in) :: DONTMASK ! If set ignore baselinequantity mask
      logical, intent(in) :: IGNORETEMPLATE
      ! Local variables
      integer :: MIF
      integer :: CHAN
      integer :: IND                    ! Combined MIF/CHAN
      integer :: i
      integer :: numProfs
      logical :: skipMe

      ! Executable code

      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ApplyBaseline', key )
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%quantityType /= l_radiance ) &
          & call Announce_Error ( key, no_Error_Code, &
          &   'Quantity to fill must be a radiance' )
        if ( baselineQuantity%template%quantityType /= l_l1bMAFBaseline ) &
          & call Announce_Error ( key, no_Error_Code, &
          &   'Quantity to fill must be a L1BMAFBaseline' )
        if ( baselineQuantity%template%signal /= quantity%template%signal .or. &
          &  baselineQuantity%template%sideband /= quantity%template%sideband ) &
          & call Announce_Error ( key, no_Error_Code, &
          &   'Quantity and baselineQuantity must have matching signal/sideband' )
      endif
      ind = 1
      numProfs = size( quantity%values ( ind, : ) )
      if ( quadrature ) then
        do mif = 1, quantity%template%noSurfs
          do chan = 1, quantity%template%noChans
            if ( .not. dontMask .and. associated(baselineQuantity%mask) ) then
              do i=1, numProfs
                skipMe = .not. dontMask .and. &
                  &  isVectorQtyMasked(baselineQuantity, chan, i, m_linAlg)
                if ( .not. skipMe )  &
                & quantity%values ( ind, i ) = sqrt ( &
                  & quantity%values ( ind, i )**2 + &
                  & baselineQuantity%values ( chan, i )**2 )
              enddo
            else
              quantity%values ( ind, : ) = sqrt ( quantity%values ( ind, : )**2 + &
                & baselineQuantity%values ( chan, : )**2 )
            endif
            ind = ind + 1
          end do
        end do
      else
        do mif = 1, quantity%template%noSurfs
          do chan = 1, quantity%template%noChans
            if ( .not. dontMask .and. associated(baselineQuantity%mask) ) then
              do i=1, numProfs
                skipMe = .not. dontMask .and. &
                  &  isVectorQtyMasked(baselineQuantity, chan, i, m_linalg)
                if ( .not. skipMe )  &
                & quantity%values ( ind, i ) = &
                  & quantity%values ( ind, i ) + &
                  & baselineQuantity%values ( chan, i )
              enddo
            else
              quantity%values ( ind, : ) = quantity%values ( ind, : ) + &
                & baselineQuantity%values ( chan, : )
            endif
            ind = ind + 1
          end do
        end do
      end if
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ApplyBaseline' )
    end subroutine ApplyBaseline

    ! ---------------------------------------------- AutoFillVector -----
    subroutine AutoFillVector ( vector )
      ! Automatically Fill items in vector we know how to
      type (Vector_T), intent(inout) :: Vector

      ! Local variables
      type (VectorValue_T), pointer :: SQ ! vector quantity
      integer :: MOL                      ! Molecule index 
      integer :: SQI                      ! Quantity index 
      character(len=32) :: str
      ! DEEBUG = .true.
      ! Executable code

      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.AutoFillVector' )
      ! Loop over its qtys
      do sqi = 1, size ( vector%quantities )
        sq => vector%quantities(sqi)
        select case(sq%template%quantityType)
        case ( l_isotoperatio )
          mol = sq%template%molecule
          if ( mol > 0 ) then
            call Get_String ( lit_indices(mol), str, strip=.true. )
            sq%values = Catalog(mol)%DefaultIsotopeRatio
            if ( DEEBUG ) then
              call outputNamedValue( 'molecule', str )
              call outputNamedValue( 'isotoperatio', Catalog(mol)%DefaultIsotopeRatio )
            endif
          endif
        case default
        end select
      enddo
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.AutoFillVector' )
    end subroutine AutoFillVector

    ! ------------------------------------- deallocateStuff ---
    subroutine deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
      real (r8), pointer, dimension(:) :: Zetab
      real (r8), pointer, dimension(:) :: Zetac
      real (r8), pointer, dimension(:) :: Zetai
      real (r8), pointer, dimension(:) :: Pb         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pc         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pi         ! p[i] in hPa
      ! call Deallocate_test ( pa, 'pa', ModuleName )
      call Deallocate_test ( pb, 'pb', ModuleName )
      call Deallocate_test ( pc, 'pc', ModuleName )
      ! call Deallocate_test ( pd, 'pd', ModuleName )
      call Deallocate_test ( pi, 'pi', ModuleName )
      ! call Deallocate_test ( Zetaa, 'Zetaa', ModuleName )
      call Deallocate_test ( Zetab, 'Zetab', ModuleName )
      call Deallocate_test ( Zetac, 'Zetac', ModuleName )
      ! call Deallocate_test ( Zetad, 'Zetad', ModuleName )
      call Deallocate_test ( Zetai, 'Zetai', ModuleName )
    end subroutine DeallocateStuff

    !=============================================== Explicit ==
    subroutine Explicit ( quantity, valuesNode, spreadFlag, force, &
      & globalUnit, channel, AzEl, options, FillValue, extraQuantity )

      ! This routine is called from MLSL2Fill to fill values from an explicit
      ! fill command line or as part of a compound Fill,
      ! Fill with height (range) specified
      
      ! Use (1): extraQuantity not present
      ! values node must be same shape as quantity or else /spread flag
      ! spreads scalar values over all (unmasked) quantity%values

      ! Use (2): extraQuantity present
      ! values node ignored
      ! sends ExtraQuantity%values into all (unmasked) quantity%values

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node
      logical, intent(in) :: SPREADFLAG   ! One instance given, spread to all
      logical, intent(in) :: FORCE        ! Fill in as many instances as will fit
      integer, intent(in) :: GLOBALUNIT   ! From parent vector
      integer, intent(in) :: CHANNEL
      logical, intent(in), optional :: AzEl ! Values are in [Mag, Az, El]; the
        ! desired quantity is components of Mag in the coordinate system to
        ! which Az and El are referenced.  So the number of values has to be
        ! a multiple of 3.
                                          ! (defaults to replacing all)
      ! The options are peculiar to this procedure, apart from verbose
      ! option           meaning
      ! ------           -------
      !   v              verbose
      !   e              replace only values in quantity == Value
      !   n              replace only values in quantity != Value
      !   a              replace only values at heights above specified height
      !   b              replace only values at heights below specified height
      !                   (defaults to replacing all)
      character (len=*), optional, intent(in) :: options ! E.g., '-v'
      real(r8), intent(in), optional :: FillValue
      type (VectorValue_T), optional :: EXTRAQUANTITY ! Instead of Value

      ! Local variables
      integer :: chan
      integer :: K                        ! Loop counter
      integer :: I,J                      ! Other indices
      logical :: MyAzEl
      real(kind(quantity%values)) :: myValue
      character (len=8) :: myOptions
      integer :: NoValues
      integer :: numChans
      integer :: surf
      integer :: TestUnit                 ! Unit to use
      integer, dimension(2) :: unitAsArray ! Unit for value given
      real (r8), pointer, dimension(:) :: VALUES
      real (r8), dimension(2) :: valueAsArray ! Value given
      logical :: Verbose
      character(len=128) :: whichInstances
      character(len=2) :: whichToReplace ! '/=' (.ne. fillValue), '==', or ' ' (always)

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.Explicit', valuesNode )
      myAzEl = .false.
      if ( present(azEl) ) myAzEl = azEl
      myOptions = ' '
      if ( present(options) ) myOptions = options
      whichInstances = ' '

      testUnit = quantity%template%unit
      if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit
      noValues = -1 ! if we will ignore valuesNode
      if ( .not. present(ExtraQuantity) ) noValues = nsons(valuesNode) - 1
      if ( Force ) noValues = min( noValues, &
        & quantity%template%instanceLen * quantity%template%noInstances )

      myValue = 0.
      if ( present(FillValue) ) myValue = FillValue

      whichToReplace = ' '
      if ( index(myOptions, 'e') > 0 ) then
        whichToReplace = '=='
      elseif ( index(myOptions, 'n') > 0 ) then
        whichToReplace = '/='
      end if
      verbose = ( index(myOptions, 'v') > 0 )
      
      if ( .not. present(extraQuantity) ) then
        ! Check the dimensions work out OK
        if ( myAzEl .and. mod(noValues,3) /= 0 ) &
            & call Announce_Error ( valuesNode, invalidExplicitFill )
        if ( spreadFlag ) then
          if ( noValues /= quantity%template%instanceLen .and. &
            & noValues /= quantity%template%noChans .and. &
            & noValues /= 1 ) &
            & call Announce_Error ( valuesNode, invalidExplicitFill )
        elseif ( .not. Force ) then
          if ( noValues /= &
            & quantity%template%instanceLen * quantity%template%noInstances ) &
            & call Announce_Error ( valuesNode, invalidExplicitFill, &
              & extraInfo = (/ &
                & quantity%template%instanceLen * quantity%template%noInstances /) )
        end if

        ! Get the values the user asked for, checking their units
        nullify ( values )
        call Allocate_test ( values, noValues, 'values', ModuleName )
        if ( .not. myAzEl ) then
          do k = 1, noValues
            call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
              & (/testUnit, PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
              & extraInfo=(/unitAsArray(1), testUnit, PHYQ_Dimensionless/) )
            values ( k ) = valueAsArray(1)
          end do
        else
          ! Convert from Mag, Az, El to 3-D projections
          do k = 1, noValues, 3
            call expr_check ( subtree(k+1,valuesNode) , unitAsArray, valueAsArray, &
              & (/testUnit, PHYQ_Dimensionless/), unitsError )
            if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
              & extraInfo=(/unitAsArray(1), testUnit, PHYQ_Dimensionless/) )
            values ( k ) = valueAsArray(1)
            ! Next two quantities have to be angles
            call expr_check ( subtree(k+2,valuesNode) , unitAsArray, valueAsArray, &
              & (/PHYQ_Angle/), unitsError )
            if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
            values (k+1) = deg2rad * valueAsArray(1)
            call expr_check ( subtree(k+3,valuesNode) , unitAsArray, valueAsArray, &
              & (/PHYQ_Angle/), unitsError )
            if ( unitsError ) call Announce_error ( valuesNode, wrongUnits, &
              & extraInfo=(/unitAsArray(1), PHYQ_Angle/) )
            values (k+2) = deg2rad * valueAsArray(1)
            values(k:k+2) = values(k) * (/ cos(values(k+1))*cos(values(k+2)), &
                                           sin(values(k+1))*cos(values(k+2)), &
                                           sin(values(k+2)) /)
          end do
        end if

        if ( verbose ) then
          call output('Explicit fill of values for ', advance='no')
          call display_string ( quantity%template%name )
          call newline
          call output(values)
          call newline
        end if
      end if

      numChans = quantity%template%instanceLen / quantity%template%noSurfs
      if ( numChans /= quantity%template%noChans ) then
        call outputNamedValue( 'noSurfs', quantity%template%noSurfs )
        call outputNamedValue( 'noChans', quantity%template%noChans )
        call outputNamedValue( 'numChans', numChans )
        call outputNamedValue( 'instanceLen', quantity%template%instanceLen )
        call announce_error ( valuesNode, no_Error_Code, &
          & 'Inconsistent template instance length' )
      endif
      ! Now loop through the quantity
      k = 0
      do i = 1, quantity%template%noInstances
        j = 0
        do surf = 1, quantity%template%noSurfs
          ! Have we specified which height to fill?
          do chan = 1, numChans
            j = j + 1
            k = k + 1
            if ( associated ( quantity%mask ) ) then
              if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) /= 0 ) cycle
            end if
            select case (whichToReplace)
            case ('/=')
              if ( quantity%values(j,i) == myValue ) cycle
            case ('==')
              if ( quantity%values(j,i) /= myValue ) cycle
            end select
            ! Have we specified which channel to fill?
            if ( channel /= 0 ) then
              if ( channel /= chan ) cycle
            endif
            if ( present(extraQuantity) ) then
              quantity%values(j,i) = extraQuantity%values(j,i)
            else
              quantity%values(j,i) = values ( mod ( k-1, noValues ) + 1 )
            endif
          end do
        end do
      end do

      if ( verbose ) then
        call output(quantity%values(1,:))
        call newline
      end if
      ! Housekeeping
      if ( .not. present(extraQuantity) ) &
        & call Deallocate_test ( values, 'values', ModuleName )
      ! No, don't do this, because sourceHeights might be our vgrid
      ! call Deallocate_test ( sourceHeights, 'sourceHeights', ModuleName )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.Explicit' )

    end subroutine Explicit

    ! ------------------------------------------- ComputeTotalPower
    subroutine ComputeTotalPower ( key, vectors )
      use MoreTree, only: GET_FIELD_ID

      ! Arguments
      integer, intent(in) :: KEY        ! Tree node
      type (Vector_T), pointer, dimension(:) :: vectors

      ! Local variables
      integer :: J                      ! Tree index
      integer :: FIELD                  ! Entry in tree
      integer :: SON                    ! Tree node
      integer :: I                      ! Loop counter
      integer :: K                      ! Loop counter
      integer :: SIGNAL                 ! the signal we're looking for

      type (Vector_T), pointer :: MEASUREMENTS ! The measurement vector to compute total power for
      type (Vector_T), pointer :: TOTALPOWERVECTOR ! The total power vector to fill
      type (Vector_T), pointer :: WEIGHTSVECTOR ! The vector of weights to use

      type (VectorValue_T), pointer :: RADIANCES ! The radiances for this band
      type (VectorValue_T), pointer :: THISRESULT ! The total power quantity to fill
      type (VectorValue_T), pointer :: WEIGHTSQUANTITY ! The total power quantity to fill

      real(rv) :: TOTALWEIGHT
      
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ComputeTotalPower', key )
      do i = 2, nsons(key)
        son = subtree(i,key)
        field = get_field_id(son)
        select case ( field )
        case ( f_measurements )
          measurements => vectors(decoration(decoration(subtree(2,son))))
        case ( f_totalPowerVector )
          totalPowerVector => vectors(decoration(decoration(subtree(2,son))))
        case ( f_weightsVector )
          weightsVector => vectors(decoration(decoration(subtree(2,son))))
        end select
      end do ! i = 2, nsons(key)

      ! Loop over the quantities in the total power vector

      do i = 1, totalPowerVector%template%noQuantities
        if ( totalPowerVector%quantities(i)%template%quantityType /= l_baseline .or. &
          & .not. totalPowerVector%quantities(i)%template%minorFrame ) then
          call Announce_Error ( key, no_error_code, 'Total power quantity must be minor frame baseline' )
        end if
        totalWeight = 0.0_rv
        thisResult => totalPowerVector%quantities(i)
        thisResult%values = 0.0_rv
        ! Now go through all the bands in the measurement vector that are in this radiometer
        do j = 1, measurements%template%noQuantities
          if ( measurements%quantities(j)%template%quantityType /= l_radiance ) cycle
          if ( measurements%quantities(j)%template%radiometer /= thisResult%template%radiometer ) cycle
          radiances => measurements%quantities(j)
          signal = measurements%quantities(j)%template%signal
          ! Now look for the weights for this band
          weightsQuantity => GetVectorQuantityByType ( weightsVector, quantityType=l_totalPowerWeight, &
            & signal=signal )
          do k = 1, radiances%template%noChans
            thisResult%values(:,:) = thisResult%values(:,:) + &
              & weightsQuantity%values(k,1) * &
              &   radiances%values(k:radiances%template%instanceLen:radiances%template%noChans,:)
            totalWeight = totalWeight + weightsQuantity%values(k,1)
          end do                        ! End loop over channels
        end do                          ! End loop over bands

!debug
!call dump(thisResult%values, 'beforedivide')


        if ( totalWeight == 0 ) then
          thisResult%values = 0         ! Don't divide by zero
        else
          thisResult%values = thisResult%values / totalWeight
        end if
      end do                            ! End loop over (effectively) radiometers
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ComputeTotalPower' )
    end subroutine ComputeTotalPower

    ! ------------------------------------------- ExtractSingleChannel ---
    subroutine ExtractSingleChannel ( key, quantity, sourceQuantity, &
      & channel, ignoreTemplate )
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Source quantity for radiances
      integer, intent(in) :: CHANNEL    ! Channel number
      logical, intent(in)              :: ignoreTemplate
      ! Local variables
      integer :: CHANIND                ! Channel index
      integer :: MIF                    ! Minor frame index
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ExtractSingleChannel', key )
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%quantityType /= l_singleChannelRadiance ) &
          & call announce_error ( key, no_Error_Code, 'Quantity to fill must be of type singleChannelRadiance' )
        if ( all ( sourceQuantity%template%quantityType /= (/ l_cloudInducedRadiance, l_radiance /) ) ) &
          & call announce_error ( key, no_Error_Code, 'source quantity for fill must be of type [cloudInduced]radiance' )
        if ( quantity%template%signal /= sourceQuantity%template%signal .or. &
          & quantity%template%sideband /= sourceQuantity%template%sideband ) &
          & call announce_error ( key, no_Error_Code, 'quantity/sourceQuantity must be same signal/sideband' )
        if ( .not. sourceQuantity%template%regular ) &
          & call Announce_Error ( key, no_Error_Code, 'source quantity must be regular' )
      endif
      chanInd = channel - GetFirstChannel ( quantity%template%signal ) + 1
      do mif = 1, quantity%template%noSurfs
        quantity%values ( mif, : ) = &
          & sourceQuantity%values ( chanInd + ( mif - 1 ) * sourceQuantity%template%noChans, : )
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ExtractSingleChannel' )
    end subroutine ExtractSingleChannel

    ! ------------------------------------------- ChiSqChan ---
    subroutine ChiSqChan ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared
      ! broken out according to channels
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in) ::      multiplier
      logical, intent(in)           ::       IGNORETEMPLATE

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             C           ! Channel loop counter
      integer ::                             S           ! Surface loop counter
      integer ::                             I           ! Instances
      integer ::                             QINDEX
      integer ::                             NOCHANS
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe
      real                             ::    a, b

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ChiSqChan', key )

      ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( &
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if ( ignoreTemplate ) then
        ! Anything goes
      elseif ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqChan/), majorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 channelwise'  )
        if ( DEEBUG ) then
          call output('major frame? ', advance = 'no')
          call output(qty%template%majorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqChan, advance = 'yes')
        end if
        go to 9
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 channelwise'  )
        go to 9
      else if ( any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing error filling chi^2 channelwise'  )
        go to 9
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) go to 9

      call allocate_test(values, measQty%template%noSurfs, &
        & 'chi^2 unsummed', ModuleName)
      noChans = qty%template%noChans
      do i=useFirstInstance, useLastInstance
        do c=1, noChans
          N = 0
          values = 0.0
          do s=1, measQty%template%noSurfs
            qIndex = c + (s-1)*nochans
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, qIndex, i, m_linalg) .or. &
            &   isVectorQtyMasked(modelQty, qIndex, i, m_linalg) .or. &
            &   isVectorQtyMasked(noiseQty, qIndex, i, m_linalg) ) &
            & .or. (ignoreNegative .and. noiseQty%values(qIndex, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(qIndex, i) == 0.0 )
            if ( .not. skipMe ) then
              values(s) = ( &
              & (a*measQty%values(qIndex, i) - b*modelQty%values(qIndex, i)) &
              & / &
              & noiseQty%values(qIndex, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(c, i) = sum(values) / N
          else
            qty%values(c, i) = 0.
          end if
        end do
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ChiSqChan' )
    end subroutine ChiSqChan

    ! ------------------------------------------- ChiSqMMaf ---
    subroutine ChiSqMMaf ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared
      ! broken out according to major frames
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in) ::      multiplier
      logical, intent(in)           ::       IGNORETEMPLATE

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             I           ! Instances
      integer ::                             ROW         ! Running 1st coord
      integer ::                             INSTANCELEN ! Num of rows
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe

      real                             ::    a, b

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ChiSqMMaf', key )

      ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( &
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if ( ignoreTemplate ) then
        ! Anything goes
      elseif ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqMMaf/), majorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 MMAFwise'  )
        if ( DEEBUG ) then
          call output('major frame? ', advance = 'no')
          call output(qty%template%majorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqMMaf, advance = 'yes')
        end if
        go to 9
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMAFwise'  )
        go to 9
      else if ( any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing noise filling chi^2 MMAFwise'  )
        go to 9
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) go to 9

      instanceLen = measQty%template%noChans * measQty%template%noSurfs
      call allocate_test(values, instanceLen, &
        & 'chi^2 unsummed', ModuleName)
      do i=useFirstInstance, useLastInstance
        if ( .not. (.not. dontMask .or. ignoreNegative .or. ignoreZero ) ) then
            values = ( &
            & (measQty%values(:, i) - modelQty%values(:, i)) &
            & / &
            & noiseQty%values(:, i) &
            &  ) ** 2
            qty%values(1, i) = sum(values) / instanceLen
        else
          N = 0
          values = 0.0
          do row = 1, instanceLen
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, row, i, m_linalg) .or. &
            &   isVectorQtyMasked(modelQty, row, i, m_linalg) .or. &
            &   isVectorQtyMasked(noiseQty, row, i, m_linalg) ) &
            & .or. (ignoreNegative .and. noiseQty%values(row, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(row, i) == 0.0 )
            if ( .not. skipMe ) then
              values(row) = ( &
              & (a*measQty%values(row, i) - b*modelQty%values(row, i)) &
              & / &
              & noiseQty%values(row, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(1, i) = sum(values) / N
          else
            qty%values(1, i) = 0.
          end if
        end if
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ChiSqMMaf', key )
    end subroutine ChiSqMMaf

    ! ------------------------------------------- ChiSqMMif ---
    subroutine ChiSqMMif ( key, qty, measQty, modelQty, noiseQty, &
    & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier, &
    & firstInstance, lastInstance )
      ! A special fill of chi squared
      ! broken out according to Mifs
      ! Formal arguments
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) ::    modelQty
      type (VectorValue_T), intent(in) ::    measQty
      type (VectorValue_T), intent(in) ::    noiseQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       ignoreZero  ! Ignore 0 values of noiseQty
      logical, intent(in)           ::       ignoreNegative  ! Ignore <0 values of noiseQty
      real, dimension(:), intent(in) ::      multiplier
      logical, intent(in)           ::       IGNORETEMPLATE

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the qty

      ! Local variables
      real(r8), dimension(:), pointer  ::    VALUES => NULL()
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             C           ! Channel loop counter
      integer ::                             S           ! Surface loop counter
      integer ::                             I           ! Instances
      integer ::                             QINDEX
      integer ::                             NOMIFS
      integer ::                             N           ! Num. of summed values
      logical ::                             skipMe

      real                             ::    a, b

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ChiSqMMif', key )

      ! Either multiplier = [a, b] or multiplier = 1/a if a=b are possible
      if ( &
      & multiplier(1) == UNDEFINED_VALUE .and. multiplier(2) == UNDEFINED_VALUE &
      & ) then
        a = 1.
        b = 1.
      else if ( multiplier(2) == UNDEFINED_VALUE ) then
        a = 0.
        if ( multiplier(1) /= 0.0 ) a = 1.0/multiplier(1)
        b = a
      else
        a = multiplier(1)
        b = multiplier(2)
      end if

      ! First check that things are OK.
      if ( ignoreTemplate ) then
        ! Anything goes
      elseif ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_chiSqMMif/), minorFrame=.true.) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 MMIFwise'  )
        if ( DEEBUG ) then
          call output('minor frame? ', advance = 'no')
          call output(qty%template%minorFrame, advance = 'no')
          call output('   quantity type ', advance = 'no')
          call output(qty%template%quantityType, advance = 'no')
          call output('   compared with ', advance = 'no')
          call output(l_chiSqMMif, advance = 'yes')
        end if
        go to 9
      else if ( .not. FillableChiSq ( qty, measQty, modelQty, noiseQty ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Incompatibility among vector quantities filling chi^2 MMIFwise'  )
        go to 9
      else if ( any ( noiseQty%values == 0.0) .and. &
        & .not. (ignoreZero .or. .not. dontMask) ) then
        call Announce_error ( key, No_Error_code, &
        & 'A vanishing noise filling chi^2 MMIFwise'  )
        go to 9
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) go to 9

      call allocate_test(values, measQty%template%noChans, &
        & 'chi^2 unsummed', ModuleName)
      noMIFs = measQty%template%noSurfs
      do i=useFirstInstance, useLastInstance
        do s=1, noMIFs
          N = 0
          values = 0.0
          do c=1, measQty%template%noChans
            qIndex = c + (s-1)*measQty%template%noChans
            skipMe = &
            & .not. dontMask .and. ( &
            &   isVectorQtyMasked(measQty, qIndex, i, m_linalg) .or. &
            &   isVectorQtyMasked(modelQty, qIndex, i, m_linalg) .or. &
            &   isVectorQtyMasked(noiseQty, qIndex, i, m_linalg) ) &
            & .or. (ignoreNegative .and. noiseQty%values(qIndex, i) < 0.0 ) &
            & .or. (ignoreZero .and. noiseQty%values(qIndex, i) == 0.0 )
            if ( .not. skipMe ) then
              values(c) = ( &
              & (a*measQty%values(qIndex, i) - b*modelQty%values(qIndex, i)) &
              & / &
              & noiseQty%values(qIndex, i) &
              &  ) ** 2
              N = N + 1
            end if
          end do
          if ( N > 0 ) then
            qty%values(s, i) = sum(values) / N
          else
            qty%values(s, i) = 0.
          end if
        end do
      end do
      call deallocate_test(values, &
        & 'chi^2 unsummed', ModuleName)
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ChiSqMMif' )
    end subroutine ChiSqMMif

    ! ------------------------------------------- ChiSqRatio ---
    subroutine ChiSqRatio ( key, qty, normQty, minNormQty, flagQty, &
    & dontMask, ignoreTemplate, firstInstance, lastInstance )
      ! A special fill of the ratio
      !  chi squared Norm
      ! ----------------     [iter_n, *]
      ! chi squared Min Norm
      ! where iter_n is the final iteration number
      
      ! Note the following tricks:
      ! The number of surfaces is the maximum allowed number of iterations
      ! The actual number of iterations will be less than this
      
      ! Depending on UNIFORMCHISQRATIO
      ! TRUE    all values will equal the ratio of tghe last ieration
      ! FALSE   nth value will be ratio for nth iteration, up to last one
      !           and all zero thereafter
      ! After the last iteration, all "surfaces" above this are zero-filled
      
      ! The number of instances will be the number of chunks
      ! (yes, an unfortunate fact)
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(inout) :: normQty
      type (VectorValue_T), intent(inout) :: minNormQty
      type (VectorValue_T), intent(inout) :: flagQty
      logical, intent(in)           ::       dontMask    ! Use even masked values
      logical, intent(in)           ::       IGNORETEMPLATE

      integer, intent(in), optional ::       firstInstance, lastInstance
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in qty

      ! Local variables
      integer ::                             UseFirstInstance, UseLastInstance, &
      &                                      NoOutputInstances
      integer ::                             I           ! Instances
      integer ::                             ITER        ! Instances
      integer ::                             QINDEX
      integer ::                             NOCHANS
      logical ::                             skipMe
      logical, parameter ::                  FakeData = .false.

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ChiSqRatio', key )

      ! First check that things are OK.
      if ( ignoreTemplate ) then
        ! Anything goes
      elseif ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_dnwt_chiSqRatio/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill wrong quantity with chi^2 ratio'  )
        go to 9
      elseif ( .not. ValidateVectorQuantity ( normqty, &
        & quantityType=(/l_dnwt_chiSqNorm/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong norm quantity with chi^2 ratio'  )
        go to 9
      elseif ( .not. ValidateVectorQuantity ( minnormqty, &
        & quantityType=(/l_dnwt_chiSqMinNorm/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong min norm quantity with chi^2 ratio'  )
        go to 9
      elseif ( .not. ValidateVectorQuantity ( flagqty, &
        & quantityType=(/l_dnwt_flag/) ) ) then
        call Announce_error ( key, No_Error_code, &
        & 'Attempting to fill using wrong flag quantity with chi^2 ratio'  )
        go to 9
      end if

      ! Work out what to do with the first and last Instance information

      if ( PRESENT(firstInstance) ) then
        useFirstInstance = firstInstance
      else
        useFirstInstance = 1
      end if

      if ( PRESENT(lastInstance) ) then
        useLastInstance = lastInstance
      else
        useLastInstance = qty%template%noInstances
      end if
      noOutputInstances = useLastInstance-useFirstInstance+1
      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) go to 9

      noChans = qty%template%noChans
      do i=useFirstInstance, useLastInstance
        if ( FakeData ) then
          ! Let's just fake up some data here
          flagQty%values(:,i) = 0._rv
          normQty%values(:,i) = 0._rv
          minNormQty%values(:,i) = 0._rv
          ! Say we converged at iteration number 5
          qIndex = 5
          flagQty%values(:qIndex,i) = -2._rv
          minNormQty%values(:qIndex,i) = -2._rv
          do iter = 1, qIndex
            minNormQty%values(iter, i) = 1.1
            normQty%values(iter, i) = 1.1 + 0.05*(qIndex+1-iter)
          enddo
        endif
        ! Now find the iteration number
        qIndex = findLast( flagQty%values(:,i)            /= 0._rv  .and. &
          &                minNormQty%values(:, i) /= 0._rv )
        if ( qIndex == 0 .or. qIndex >= qty%template%noSurfs ) cycle
        skipMe = &
          & .not. dontMask .and. ( &
          &   isVectorQtyMasked(normQty, qIndex, i, m_linalg) .or. &
          &   isVectorQtyMasked(minNormQty, qIndex, i, m_linalg) .or. &
          &   minNormQty%values(qIndex, i) == 0. &
          & )
        if ( skipMe ) then
          if ( minNormQty%values(qIndex, i) == 0. ) qty%values(:,i) = 999._rv
        elseif ( UNIFORMCHISQRATIO .or. &
          & size(qty%values) /= size(normQty%values) .or. &
          & size(qty%values) /= size(minNormQty%values) ) then
          qty%values(:,i) = &
            & normQty%values(qIndex, i) / minNormQty%values(qIndex, i)
        else
          qty%values(:,i) = 0._rv
          qty%values(1:qIndex,i) = &
            & normQty%values(1:qIndex, i) / minNormQty%values(1:qIndex, i)
        endif
      end do
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ChiSqRatio' )
    end subroutine ChiSqRatio

    ! ------------------------------------------- ColAbundance ---
    subroutine ColAbundance ( key, qty, bndPressQty, vmrQty, &
      & colmAbUnits, ignoreTemplate, &
      & firstInstance, lastInstance )
      ! A special fill according to W.R.Read's idl code
      ! Similar to his hand-written notes, but with a small correction

      ! Assumptions:
      ! (See above)
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: BNDPRESSQTY
      type (VectorValue_T), intent(in) :: VMRQTY
      integer, intent(in) :: colmAbUnits
      logical, intent(in)           :: IGNORETEMPLATE
      integer, intent(in), optional :: FIRSTINSTANCE
      integer, intent(in), optional :: LASTINSTANCE
      ! The last two are set if only part (e.g. overlap regions) of the quantity
      ! is to be stored in the column data.

      ! Local parameters
      real(r8), parameter :: INVERMGPPMV = 0.789 ! in DU / (ppmv hPa)
      real(r8), parameter :: INVERMGDU = INVERMGPPMV*1.d6 ! in DU / (vmr hPa)
      real(r8), parameter :: INVERMGMOLCM2 = INVERMGDU*2.687d16 ! in mol / cm^2

      ! Local variables
      logical :: printMe
      integer :: SURFACE
      integer :: INSTANCE
      integer :: FIRSTSURFACEABOVE
      integer :: FIRSTSURFACEBELOW
      integer :: USEFIRSTINSTANCE
      integer :: USELASTINSTANCE
      integer :: N
      integer :: NOOUTPUTINSTANCES
      real (r8) :: THISBNDPRESS
      real (r8) :: zeta            ! -log(THISBNDPRESS)
      real (r8) :: zetaTmp
      real (r8) :: COLUMNSUM
      real (r8) :: TRAPEZOIDSUM
      real (r8) :: DELTAZETA      ! Zetai[s+1] - Zetai[s]
      real (r8) :: INVERMG
      real (r8)                        :: Zetaa
      real (r8), pointer, dimension(:) :: Zetab
      real (r8), pointer, dimension(:) :: Zetac
      real (r8)                        :: Zetad
      real (r8), pointer, dimension(:) :: Zetai
      real (r8)                        :: Pa         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pb         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pc         ! p[i] in hPa
      real (r8)                        :: Pd         ! p[i] in hPa
      real (r8), pointer, dimension(:) :: Pi         ! p[i] in hPa

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ColAbundance', key )
      ! First check that things are OK.
      if ( ignoreTemplate ) then
        ! Anything goes
      elseif ( (qty%template%quantityType /= l_columnAbundance) .or. &
        &  (bndPressQty%template%quantityType /= l_boundaryPressure) .or. &
        &  (vmrQty%template%quantityType /= l_vmr) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Wrong quantity type found while filling column abundance'  )
        go to 9
      else if ( qty%template%molecule /= vmrQty%template%molecule ) then
        call Announce_error ( key, No_Error_code, &
          & 'Attempt to fill column abundance with different molecule'  )
        go to 9
      else if ( .not. ( DoHgridsMatch( qty, vmrQty ) .and. &
        & DoHgridsMatch( qty, bndPressQty ) ) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Attempt to fill column abundance with different HGrids'  )
        go to 9
      else if ( .not. any(vmrQty%template%verticalCoordinate == &
        & (/l_zeta/)) ) then
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but vmr not on zeta surfs.'  )
        go to 9
      else if ( vmrQty%template%noSurfs < 2 ) then
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but too few vmr surfaces'  )
        go to 9
      end if
      select case (colmAbUnits)
      case (l_DobsonUnits)
        INVERMG = INVERMGDU
      case (l_DU)
        INVERMG = INVERMGDU
      case (l_molcm2)
        INVERMG = INVERMGMOLCM2
      case default
        call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, but wrong units.'  )
      end select
      ! Work out what to do with the first and last Instance information
      useFirstInstance = 1
      useLastInstance = qty%template%noInstances
      if ( present ( firstInstance ) ) useFirstInstance = firstInstance
      if ( present ( lastInstance ) ) useLastInstance = lastInstance
      noOutputInstances = useLastInstance-useFirstInstance+1

      ! If we've not been asked to output anything then don't carry on
      if ( noOutputInstances < 1 ) go to 9

      nullify ( Zetab, Zetac, Zetai, pb, pc, pi )
      call allocate_test ( Zetab, vmrQty%template%noSurfs, 'Zetab', ModuleName )
      call allocate_test ( Zetac, vmrQty%template%noSurfs, 'Zetac', ModuleName )
      call allocate_test ( Zetai, vmrQty%template%noSurfs, 'Zetai', ModuleName )
      call allocate_test ( pb, vmrQty%template%noSurfs, 'pb', ModuleName )
      call allocate_test ( pc, vmrQty%template%noSurfs, 'pc', ModuleName )
      call allocate_test ( pi, vmrQty%template%noSurfs, 'pi', ModuleName )
      Zetai = vmrQty%template%surfs(:,1)
      pi = 10.0 ** ( -Zetai  )
      N = vmrQty%template%noSurfs
      do instance = useFirstInstance, useLastInstance
        printMe = ( instance == useFirstInstance) .and. &
         & (switchDetail(switches, 'column') > -1 )
        if ( printMe ) print *, 'switches: ', trim(switches)
        if ( printMe ) &
          & print *, 'switchDetail(switches, column) ', switchDetail(switches, 'column')
        ! Find 1st surface at or above tropopause
        ! (i.e., at a pressure equal to or less than boundaryPressure)
        ! This next check should be unnecessary--the HGrids were already matched
        if ( instance > size(bndPressQty%values, 2) ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Cant fill column--instance outside b.pres. range' )
          call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
          go to 9
        end if
        thisBndPress = bndPressQty%values(1,instance)
        ! In case where WMO algorithm failed, use bottom of basis
        if ( thisBndPress <= 0.0 ) &
          & thisBndPress = 10.0 ** ( - vmrQty%template%surfs(1,1) )
        if ( thisBndPress <= 0._r8 ) then
          call Announce_error ( key, No_Error_code, &
          & 'Fill column abundance, illegal bound. pr. at this instance' )
        end if
        zeta = -log10 ( thisBndPress )
        Zetaa = max(Zetai(N), zeta)
        do surface=1, N
          Zetab(surface) = Zetai(min(surface+1,N))
        end do
        do surface=1, N
          Zetab(surface) = min(max(zeta, Zetai(surface)), Zetab(surface))
          zetaTmp = zeta
          if ( surface > 1 ) zetaTmp = max(zeta, Zetai(surface-1))
          Zetac(surface) = min(zetaTmp, Zetai(surface))
        end do
        Zetad = min(Zetai(1), zeta)
        Pa = 10. ** (-Zetaa)
        Pb = 10. ** (-Zetab)
        Pc = 10. ** (-Zetac)
        Pd = 10. ** (-Zetad)
        if ( printMe ) then
          call output ( 'zeta: ', advance='no')
          call output ( zeta, advance='yes')
          call output ( 'zetaA: ', advance='no')
          call output ( zetaA, advance='yes')
          call output ( 'zetaD: ', advance='no')
          call output ( zetaD, advance='yes')
          call output ( 'PA: ', advance='no')
          call output ( PA, advance='yes')
          call output ( 'PD: ', advance='no')
          call output ( PD, advance='yes')
          call output ( 'zetaI   zetaB     zetaC      Pi   Pb   Pc', advance='yes')
          do surface=1, N
            call output( (/ zetai(surface), zetab(surface), zetac(surface), &
             & Pi(surface), Pb(surface), Pc(surface)/) , advance='yes')
          end do
        end if
        ! Find 1st surface immediately above tropopause
        firstSurfaceAbove = FindFirst (Pi < thisBndPress)
        if ( firstSurfaceAbove < 1 ) then
          call output ( 'tropopause: ', advance='no')
          call output ( thisBndPress, advance='yes')
          call output ( 'p0 ', advance='no')
          call output ( Pi(1), advance='yes')
          call output ( 'pTop ', advance='no')
          call output ( Pi(N), advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Filling column, but tropopause outside pres. surfaces' )
          firstSurfaceBelow = 1
        else if ( firstSurfaceAbove == 1 ) then
          ! Nothing special
          firstSurfaceBelow = 1
        else if ( firstSurfaceAbove > N ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Cant column, tropopause above top surface' )
          call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
          go to 9
        else  ! Nothing special
          firstSurfaceBelow = firstSurfaceAbove - 1
        end if
        if ( printMe ) then
          print *, 'thisBndPress: ', thisBndPress
          print *, 'firstSurfaceAbove: ', firstSurfaceAbove
          print *, 'firstSurfaceBelow: ', firstSurfaceBelow
        end if
        ! Do summation
        columnSum = vmrQty%values(N, instance) * Pa  ! Initialize sum
        if ( printMe ) then
          print *, 'columnSum: ', columnSum
          trapezoidSum = 0.
          do surface=firstSurfaceAbove, N-1
            trapezoidsum = trapezoidsum + 0.5 * ( &
              & vmrQty%values(surface, instance) + vmrQty%values(surface+1, instance) &
              & ) * ( Pi(surface+1)-Pi(surface))
          end do
        end if
        ! Loop over surfaces from 1 to uppermost-1
        do surface = 1, N-1
          deltaZeta = Zetai(surface+1) - Zetai(surface)
          if ( deltaZeta == 0._r8 ) cycle
          columnSum = columnSum + &
            & vmrQty%values(surface, instance) / deltaZeta * &
            & ( &
            & Pb(surface)*(Zetai(surface+1)-Zetab(surface)) &
            & + &
            & (Pi(surface+1) - Pb(surface))/ln10 &
            & )
        if ( printMe ) print *, 'columnSum: ', columnSum, surface
        end do
        ! Loop over surfaces from 2 to uppermost
        do surface = 2, N
          deltaZeta = Zetai(surface) - Zetai(surface-1)
          if ( deltaZeta == 0._r8 ) cycle
          columnSum = columnSum + &
            & vmrQty%values(surface, instance) / deltaZeta * &
            & ( &
            & Pc(surface)*(Zetac(surface)-Zetai(surface)) &
            & + &
            & (Pc(surface) - Pi(surface))*(1./ln10 + deltaZeta) &
            & )
          ! columnSum = columnSum + &
          ! & vmrQty%values(surface, instance) / deltaZeta * &
          ! & ( &
          ! & Pc(surface)*(Zetac(surface)-Zetai(surface)) &
          ! & + &
          ! & (Pc(surface) - Pi(surface))/ln10 &
          ! & )
        if ( printMe ) print *, 'columnSum: ', columnSum, surface
        end do
        columnSum = columnSum + vmrQty%values(1, instance) * (Pd - Pi(1))
        if ( printMe ) print *, 'columnSum: ', columnSum
        if ( printMe ) print *, 'trapezoid: ', -trapezoidSum
        qty%values ( 1, instance ) = InverMg * columnSum
      end do
      call deallocateStuff(Zetab, Zetac, Zetai, Pb, Pc, Pi)
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ColAbundance' )
    end subroutine ColAbundance

    ! -------------------------------------------- ConvergenceFromChisq --------
    subroutine ConvergenceFromChisq ( key, quantity, sourceQuantity, &
      & scale, ignoreTemplate )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in) :: SOURCEQUANTITY ! dnwt_ChisqRatio quantity on which it's based
      real(r8), intent(in) :: SCALE     ! A scale factor
      logical, intent(in)           ::       IGNORETEMPLATE
      ! Local variables
      integer ::                             QINDEX
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ConvergenceFromChisq', key )
      ! Do some sanity checking
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%quantityType /= l_quality ) call Announce_error ( key, no_error_code, &
          & 'Convergence quantity must be quality' )
        if ( sourceQuantity%template%quantityType /= l_dnwt_chisqRatio ) call Announce_error ( &
          & key, no_error_code, 'sourceQuantity must be of type chisqRatio' )
      endif
      if ( UNIFORMCHISQRATIO ) then
        quantity%values(1,:) = scale * sourceQuantity%values(1,1)
      else
        qIndex = findLast( sourceQuantity%values(:,1) /= 0._rv )
        if ( qIndex > 0 ) &
          & quantity%values(1,:) = scale * sourceQuantity%values(qIndex,1)
      endif
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ConvergenceFromChisq' )
    end subroutine ConvergenceFromChisq

    !------------------------------------- FillCovariance ------------
    subroutine FillCovariance ( covariance, vectors, diagonal, &
      & lengthScale, fraction, invert, ignoreTemplate )
      ! This routine fills a covariance matrix from a given set of vectors
      type (Matrix_SPD_T), intent(inout) :: COVARIANCE ! The matrix to fill
      type (Vector_T), dimension(:), intent(in), target :: VECTORS ! The vector database
      integer, intent(in) :: DIAGONAL     ! Index of vector describing diagonal
      integer, intent(in) :: LENGTHSCALE  ! Index of vector describing length scale
      integer, intent(in) :: FRACTION     ! Index of vector describing fraction
      logical, intent(in) :: INVERT       ! We actually want the inverse
      logical, intent(in) :: IGNORETEMPLATE

      ! Local parameters
      real(r8), parameter :: DECADE = 16000.0 ! Number of meters per decade.

      ! Local variables
      type (VectorValue_T), pointer :: d  ! Diagonal
      type (VectorValue_T), pointer :: l  ! Length
      type (VectorValue_T), pointer :: f  ! Fraction
      type (QuantityTemplate_T), pointer :: qt ! One quantity template
      type (Vector_T) :: DMASKED ! Masked diagonal
      type (Vector_T) :: LMASKED ! Masked length scale
      integer :: B                        ! Block index
      integer :: I                        ! Instance index
      integer :: J                        ! Loop index
      integer :: K                        ! Loop index
      integer :: N                        ! Size of matrix block
      integer :: Q                        ! Quantity index
      type (MatrixElement_t), pointer :: M ! The matrix being filled
      real (r8), dimension(:), pointer :: SURFS ! The vertical coordinate
      real (r8) :: distance               ! Distance between two points
      real (r8) :: thisLength             ! Geometric mean length scale
      real (r8) :: meanDiag               ! Geometric mean diagonal value
      real (r8) :: thisFraction           ! Geometric mean diagonal value
      logical, dimension(:), pointer :: condition ! Condition
      logical :: ANYOFFDIAG             ! Flag to indicate presence of off diagonal elements

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FillCovariance' )

      ! Apply mask to diagonal
      nullify ( condition )
      call CopyVector ( Dmasked, vectors(diagonal), clone=.true., &
        & vectorNameText='_Dmasked' )
      call ClearUnderMask ( Dmasked )

      if ( lengthScale == 0 ) then
        call updateDiagonal ( covariance, vectors(diagonal), square=.true.,&
          & invert=invert, forgiveZeros=.true. )
      else
        ! Do a more complex fill, either we're doing non-diagonal, or there might
        ! be zeros to 'invert'

        ! Setup some stuff
        if ( lengthScale /= 0 ) then
          call CopyVector ( Lmasked, vectors(lengthScale), clone=.true., &
            & vectorNameText='_Lmasked' )
          call ClearUnderMask ( Lmasked )
        end if

        ! Check the validity of the supplied vectors
        if ( .not. ignoretemplate ) then
          if ( covariance%m%row%vec%template%name /= &
            & vectors(diagonal)%template%name ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "Diagonal and covariance not compatible in fillCovariance" )
          if ( covariance%m%row%vec%template%name /= &
            & dMasked%template%name ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "Copied diagonal and covariance not compatible in fillCovariance" )
          if ( lengthScale /= 0 ) then    ! Check length if supplied
            if ( covariance%m%row%vec%template%name /= &
              & lMasked%template%name ) call MLSMessage ( MLSMSG_Error, &
              & ModuleName, "lengthScale and covariance not compatible in fillCovariance" )
            if ( lMasked%globalUnit /= phyq_length ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "length vector does not have dimensions of length" )
          else
            thisLength = 0.0
          end if
          if ( fraction /= 0 ) then       ! Check fraction if supplied
            if ( covariance%m%row%vec%template%name /= &
              & vectors(fraction)%template%name ) call MLSMessage ( MLSMSG_Error, &
              & ModuleName, "fraction and covariance not compatible in fillCovariance" )
            if ( vectors(fraction)%globalUnit /= phyq_dimensionless ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "fraction vector is not dimensionless" )
          else
            thisFraction = 1.0
          end if
        end if

        ! Now loop over the quantities
        do q = 1, covariance%m%col%vec%template%noQuantities

          ! Setup pointers etc.
          d => dMasked%quantities(q)
          qt => d%template
          if ( lengthScale /= 0 ) l => lMasked%quantities(q)
          if ( fraction /=0 ) f => vectors(fraction)%quantities(q)
          n = qt%instanceLen
          if ( qt%coherent ) surfs => qt%surfs(:,1)
          if ( .not. qt%regular ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Unable to handle irregular quantity in Covariance" )

          ! Loop over the instances
          do i = 1, qt%noInstances
            b = FindBlock ( covariance%m%col, q, i )
            m => covariance%m%block(b,b)
            call createBlock ( m, n, n, m_full, forWhom=moduleName )
            anyOffDiag = .false.
            if ( .not. qt%coherent ) surfs => qt%surfs(:,i)

            ! Clear the working matrix and load the diagonal
            m%values = 0.0_rm
            do j = 1, n
              m%values(j,j) = d%values(j,i) ** 2.0
            end do

            ! Now if appropriate add off diagonal terms.
            if ( any( qt%verticalCoordinate == (/ l_height, l_pressure, l_zeta /) ) ) then
              ! Loop over off diagonal terms
              do j = 1, n
                do k = 1, j-1
                  meanDiag = sqrt ( m%values(j,j) * m%values(k,k) )
                  if ( lengthScale /= 0 ) &
                    & thisLength = sqrt ( l%values(j,i) * l%values(k,i) )
                  if ( fraction /= 0 ) thisFraction = f%values(j,i)
                  select case (qt%verticalCoordinate)
                  case ( l_height )
                    distance = abs ( surfs ( (j-1)/qt%noChans + 1 ) - &
                      & surfs ( (k-1)/qt%noChans + 1) )
                  case ( l_zeta )
                    distance = abs ( surfs ( (j-1)/qt%noChans + 1 ) - &
                      & surfs ( (k-1)/qt%noChans + 1 ) ) * decade
                  case ( l_pressure )
                    distance = abs ( -log10 ( surfs( (j-1)/qt%noChans + 1) ) + &
                      &               log10 ( surfs( (k-1)/qt%noChans + 1) ) ) / decade
                  end select
                  if ( thisLength > 0.0 .and. thisFraction > 0.0 ) then
                    m%values(j,k) = meanDiag*thisFraction*exp(-distance/thisLength)
                    anyOffDiag = .true.
                  end if
                end do                    ! Loop over k (in M)
              end do                      ! Loop over j (in M)
            end if                        ! An appropriate vertical coordinate

            ! Now we may need to invert this, if so we need to be clever.
            if ( invert ) then
              call Allocate_test ( condition, n, 'condition', ModuleName )
              condition = d%values(:,i) <= 0.0_rv
              do j = 1, n
                if ( condition(j) ) m%values(j,j) = 1.0_rm
              end do
              if ( anyOffDiag ) then
                call MatrixInversion(M%values, upper=.true.)
              else
                do j = 1, n
                  m%values(j,j) = 1.0 / m%values(j,j)
                end do
              end if
              do j = 1, n
                if ( condition(j) ) m%values(j,j) = 0.0_rm
              end do
              call Deallocate_test ( condition, 'condition', ModuleName )
            end if

            call Sparsify ( M )
          end do                          ! Loop over instances
        end do                            ! Loop over quantities
      end if                              ! A non diagonal fill

      call DestroyVectorInfo ( DMasked )
      call DestroyVectorInfo ( LMasked )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FillCovariance' )

    end subroutine FillCovariance

    !------------------------------------- DerivativeOfSource ------------
    ! Compute derivative of source quantity w.r.t. dim listed
    subroutine DerivativeOfSource ( DERIVATIVE, SOURCE, XQUANTITY, &
      & DIMLIST, IGNORETEMPLATE )
      use MANIPULATIONUTILS, only: MANIPULATE

      type (VectorValue_T), pointer       :: DERIVATIVE
      type (VectorValue_T), pointer       :: SOURCE
      type (VectorValue_T), pointer       :: XQUANTITY
      character(len=*), intent(in)        :: DIMLIST        ! 's', 'c', or 'i' in manipulation's shi
      logical, intent(in)                 :: IGNORETEMPLATE

      ! Local variables
      type(vectorValue_T), pointer :: BOGUS         ! When we have no "b"
      integer                      :: NOKEY
      type(vectorValue_T), pointer :: LOWER         ! For storing (dy/dx)_l
      type(vectorValue_T), target  :: LOWERQUANTITY ! For storing (dy/dx)_l
      type(vectorValue_T), pointer :: UPPER         ! For storing (dy/dx)_u
      type(vectorValue_T), target  :: UPPERQUANTITY ! For storing (dy/dx)_u
      ! Executable
      nullify( bogus, upper, lower )
      nokey = 0
      ! First: check that derivative and source are compatible
      if ( .not. ignoretemplate ) then
        if ( source%template%name /= &
          & derivative%template%name ) call MLSMessage ( MLSMSG_Error, &
          & ModuleName, "Derivative and source not compatible in fill derivative" )
      endif
      call CloneVectorQuantity ( upperQuantity, source )
      call CloneVectorQuantity ( lowerQuantity, source )
      lower => lowerQuantity
      upper => upperQuantity
      ! Second: compute the "upper derivative"
      ! (dy/dx)_u = ( y[x+dx] - y[x] ) / dx
      call Manipulate ( upper, source, xQuantity, 0._rv, &
        & "(shift(a) - a)/(shift(b) - b)", &
        & .false., dimList )
      ! Third: compute the "lower derivative"
      ! (dy/dx)_l = ( y[x] - y[x-dx] ) / dx
      call Manipulate ( lower, source, xQuantity, 0._rv, &
        & "(a - slip(a))/(b - slip(b))", &
        & .false., dimList )
      ! Last: average the two to arrive at the central difference
      ! (dy/dx) = (1/2) * ( (dy/dx)_u + (dy/dx)_l )
      call Manipulate ( derivative, upper, lower, 0.5_rv, &
        & "c*(a + b)", &
        & .false., dimList )
      ! Housekeeping
      call destroyVectorQuantityValue ( upperQuantity, &
        & destroyMask=.true., destroyTemplate=.false. )
      call destroyVectorQuantityValue ( lowerQuantity, &
        & destroyMask=.true., destroyTemplate=.false. )
    end subroutine DerivativeOfSource

    ! ------------------------------------- FoldedRadiance ---
    subroutine FoldedRadiance ( radiance, lsb, usb, &
      & lsbFraction, usbFraction, key )
      type (VectorValue_T), intent(inout) :: RADIANCE
      type (VectorValue_T), pointer :: USB
      type (VectorValue_T), pointer :: LSB
      type (VectorValue_T), pointer :: USBFRACTION
      type (VectorValue_T), pointer :: LSBFRACTION
      integer, intent(in) :: KEY

      ! Local variables
      integer :: C                        ! Channel loop inductor
      integer :: I                        ! Array index
      integer :: MIF                      ! Minor frame loop inductor

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FoldedRadiance', key )
      ! First some sanity checks
      if ( .not. ValidateVectorQuantity ( radiance, quantityType=(/l_radiance/), &
        & sideband=(/0/), minorFrame=.true. )) &
        & call Announce_Error ( key, no_error_code, 'Inappropriate radiance quantity to fill' )
      if ( ( associated ( lsb ) .neqv. associated ( lsbFraction ) ) .or. &
        &  ( associated ( usb ) .neqv. associated ( usbFraction ) ) ) then
        call Announce_Error ( key, no_error_code, 'Must supply sidebands and fractions together' )
        go to 9
      end if
      if ( associated ( lsb ) ) then
        if ( .not. ValidateVectorQuantity ( lsb, quantityType=(/l_radiance/), &
          & sideband=(/-1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsb radiance quantity for fill' )
        if ( .not. ValidateVectorQuantity ( lsbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/radiance%template%signal/), sideband=(/-1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsbFraction quantity for fill' )
      end if
      if ( associated ( usb ) ) then
        if ( .not. ValidateVectorQuantity ( usb, quantityType=(/l_radiance/), &
          & sideband=(/1/), signal=(/radiance%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usb radiance quantity for fill' )
        if ( .not. ValidateVectorQuantity ( usbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/radiance%template%signal/), sideband=(/1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usbFraction quantity for fill' )
      end if

      if ( .not. associated ( lsb ) .and. .not. associated ( usb ) ) &
        & call Announce_Error ( key, no_error_code, 'Must supply one or both of lsb/usb' )

      ! Now do the work
      ! Note that this is a very inefficient loop, but I don't really care as it's
      ! Never used in routine processing.
      i = 1                               ! Use i as a composit mif,channel index
      do mif = 1, radiance%template%noSurfs
        do c = 1, radiance%template%noChans
          radiance%values(i,:) = 0.0
          if ( associated ( lsb ) ) radiance%values(i,:) = radiance%values(i,:) + &
            & lsbFraction%values(c,1) * lsb%values(i,:)
          if ( associated ( usb ) ) radiance%values(i,:) = radiance%values(i,:) + &
            & usbFraction%values(c,1) * usb%values(i,:)
          i = i + 1
        end do
      end do
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FoldedRadiance' )

    end subroutine FoldedRadiance

    ! --------------------------------------------- FromAnother ---
    ! The most basic of Fill methods:
    ! At some point we must try to modify most of the other Fill methods so 
    ! that they do their unique stuff and then call this method
    ! That makes maximum reuse, a best practice, and confines future changes to
    ! just this subroutine
    ! Note that this does *NOT* copy the mask (at least for the moment)
    ! It is assumed that the original one (e.g. inherited from transfer)
    ! is still relevant.
    ! See Subset command
    
    ! Next:
    ! Check that vGrids match
    ! If they don't, and /interpolate not set, then raise an exception
    ! (no more silent inerpolations)
    subroutine FromAnother ( quantity, sourceQuantity, ptan, &
      & key, ignoreTemplate, spreadflag, interpolate, force )
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), pointer :: SOURCEQUANTITY
      type (VectorValue_T), pointer :: PTAN
      integer, intent(in) :: KEY        ! Tree node
      logical, intent(in) :: IGNORETEMPLATE ! If set throw caution to the wind
      logical, intent(in) :: SPREADFLAG ! If set spread across instances
      logical, intent(in) :: INTERPOLATE ! If set spread across summed dimension
      logical, intent(in) :: FORCE ! Copy as much as will fit
      ! Local parameters
      integer :: inst
      integer :: instanceLen
      integer :: noInstances
      integer :: surf
      ! Executable
      if ( associated ( ptan ) .and. interpolate ) then
        call FromInterpolatedQty ( quantity, sourceQuantity, &
          & key, ignoreTemplate, ptan )
      elseif ( quantity%template%name /= sourceQuantity%template%name ) then
        if ( .not. interpolate .and. .not. ignoreTemplate ) then
          call Announce_Error ( key, No_Error_Code, &
            & 'Quantity and sourceQuantity do not have the same template' )
        else
          call FromInterpolatedQty ( quantity, sourceQuantity, &
            & key, ignoreTemplate )
        end if
      elseif ( spreadFlag ) then
        ! Copy first instance into all (assuming instance lengths are the same)
        noInstances = Quantity%template%noInstances
        do inst = 1, noInstances
          ! If we have a mask and we're going to obey it then do so
          if ( associated(quantity%mask) ) then
            where ( iand ( ichar(quantity%mask(:,inst)), m_Fill ) == 0 )
              quantity%values(:,inst) = sourceQuantity%values(:,1)
            end where
          else ! Otherwise, just blindly copy
            quantity%values(:,inst) = sourceQuantity%values(:,1)
          end if
        enddo
      elseif ( force .and. &
        & sourceQuantity%template%instanceLen == Quantity%template%instanceLen ) then
        ! Copy as many instances as we can
        noInstances = min( sourceQuantity%template%noInstances, Quantity%template%noInstances )
        do inst = 1, noInstances
          ! If we have a mask and we're going to obey it then do so
          if ( associated(quantity%mask) ) then
            where ( iand ( ichar(quantity%mask(:,inst)), m_Fill ) == 0 )
              quantity%values(:,inst) = sourceQuantity%values(:,inst)
            end where
          else ! Otherwise, just blindly copy
            quantity%values(:,inst) = sourceQuantity%values(:,inst)
          end if
        enddo
      elseif ( force .and. &
        & sourceQuantity%template%noInstances == Quantity%template%noInstances ) then
        ! Copy as many instances as we can
        instanceLen = min( sourceQuantity%template%instanceLen, Quantity%template%instanceLen )
        do surf = 1, instanceLen
          ! If we have a mask and we're going to obey it then do so
          if ( associated(quantity%mask) ) then
            where ( iand ( ichar(quantity%mask(surf,:)), m_Fill ) == 0 )
              quantity%values(surf,:) = sourceQuantity%values(surf,:)
            end where
          else ! Otherwise, just blindly copy
            quantity%values(surf,:) = sourceQuantity%values(surf,:)
          end if
        enddo
      else
        ! Just a straight copy
        ! If we have a mask and we're going to obey it then do so
        if ( associated(quantity%mask) ) then
          where ( iand ( ichar(quantity%mask(:,:)), m_Fill ) == 0 )
            quantity%values(:,:) = sourceQuantity%values(:,:)
          end where
        else ! Otherwise, just blindly copy
          quantity%values = sourceQuantity%values
        end if
      end if

    end subroutine FromAnother

    ! ------------------------------------- FromSplitSideband ----
    subroutine FromSplitSideband ( quantity, sourceQuantity, &
      & lsbFraction, usbFraction, spreadFlag, usb, channel, key )

      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), intent(in) :: LSBFRACTION
      type (VectorValue_T), intent(in) :: USBFRACTION
      type (VectorValue_T), pointer :: USB
      logical, intent(in) :: SPREADFLAG   ! One instance given, spread to all
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: KEY
      ! Local variables
      integer :: MYCHANNEL              ! Possibly offset channel
      integer :: i, mif, maf
      type (Signal_T) :: signalIn, signalOut, signalRef
      real(r8), dimension(:), pointer :: freq1, freqL1, freqU1
      real(r8), dimension(:), pointer :: freq2, freqL2, freqU2
      real(r8), dimension(:), pointer :: freq, freqL, freqU
      real(r8) :: ratio1, ratio2    ! signal sideband fractions
      real(r8) :: scaledRad   ! scaled radiance according to the f^4 law

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromSplitSideband', key )
      nullify(freq, freqL, freqU, freq1, freqL1, freqU1, freq2, freqL2, freqU2)

      ! check for qualified quantity
      if ( .not. ValidateVectorQuantity ( quantity, quantityType=(/l_cloudInducedRadiance/), &
        & sideband=(/-1,1/), minorFrame=.true. )) &
        & call Announce_Error ( key, no_error_code, 'Quantity must be cloud-induced-radiances to fill' )
      if ( .not. associated(quantity%mask)) &
        & call Announce_Error ( key, no_error_code, 'Quantity must be a subset to fill' )

      signalIn = GetSignal ( sourceQuantity%template%signal )     ! sideband info gets lost
      signalOut = GetSignal ( Quantity%template%signal )          ! sideband info gets lost

      myChannel = channel

      call allocate_test ( freq1, size(signalIn%frequencies), 'frequencies', ModuleName )
      call allocate_test ( freqL1, size(signalIn%frequencies), 'LSBfrequencies', ModuleName )
      call allocate_test ( freqU1, size(signalIn%frequencies), 'USBfrequencies', ModuleName )
      call allocate_test ( freq2, size(signalOut%frequencies), 'frequencies', ModuleName )
      call allocate_test ( freqL2, size(signalOut%frequencies), 'LSBfrequencies', ModuleName )
      call allocate_test ( freqU2, size(signalOut%frequencies), 'USBfrequencies', ModuleName )

      ! find input signal frequencies
      freqL1 = signalIn%lo - signalIn%centerFrequency - &
        & signalIn%direction*signalIn%frequencies    ! lower sideband freq
      freqU1 = signalIn%lo + signalIn%centerFrequency + &
        & signalIn%direction*signalIn%frequencies     ! upper sideband freq
      if ( sourceQuantity%template%sideband == -1) freq1 = freqL1
      if ( sourceQuantity%template%sideband == 1) freq1 = freqU1
      ! find output signal frequencies
      freqL2 = signalOut%lo - signalOut%centerFrequency - &
        & signalOut%direction*signalOut%frequencies      ! lower sideband freq
      freqU2 = signalOut%lo + signalOut%centerFrequency + &
        & signalOut%direction*signalOut%frequencies    ! upper sideband freq
      if ( quantity%template%sideband == -1) freq2 = freqL2
      if ( quantity%template%sideband == 1) freq2 = freqU2

      quantity%values=0._r8
      if ( spreadFlag ) then
        ! spread a cloudy radiance to other bands according to the f^4 law
        ! The source cloudy radiance must be from a single sideband signal
        if ( .not. ValidateVectorQuantity ( sourceQuantity, &
          & quantityType=(/l_cloudInducedRadiance/), &
          & sideband=(/-1,1/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate sourceQuantity radiance for fill' )
        do i=1,size(signalOut%frequencies)
          do maf=1, quantity%template%noInstances
            do mif=1, quantity%template%noSurfs
              scaledRad = sourceQuantity%values ( MyChannel + &
                & (mif-1) * size(signalIn%frequencies), maf) * &
                & freq2(i)**4/freq1(MyChannel)**4
              if ( iand ( ichar ( Quantity%mask(i+(mif-1) * &
                & size(signalOut%frequencies), maf)), m_cloud) /= 0 ) &
                & quantity%values(i+(mif-1)*size(signalOut%frequencies), maf) = scaledRad
            end do
          end do
        end do

      else
        ! split two sideband radiances according to the f^4 law
        ! The source cloudy radiance must be from the same signal of double sideband

        if ( .not. ValidateVectorQuantity ( sourceQuantity, &
          & quantityType=(/l_cloudInducedRadiance/), &
          & sideband=(/0/), signal=(/quantity%template%signal/), minorFrame=.true. )) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate sourceQuantity radiance for fill' )
        if ( .not. ValidateVectorQuantity ( lsbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/quantity%template%signal/), sideband=(/-1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate lsbFraction quantity for fill' )
        if ( .not. ValidateVectorQuantity ( usbFraction, &
          & quantityType=(/l_limbSidebandFraction/), &
          & signal=(/quantity%template%signal/), sideband=(/1/) ) ) &
          & call Announce_Error ( key, no_error_code, 'Inappropriate usbFraction quantity for fill' )

        ! This method is only appliable to the cloud-induced radiances. One may assume
        ! that the scattering radiances are proportional to f^4. And this operation
        ! is only applied to maskbit = m_cloud.

        if ( .not. associated(usb) ) then
          ! If both sidebands are within 20GHz and have similar penetration depths
          ! in the attenuative atmosphere. We can neglect the attenuation and split cloudy
          ! radiances assuming that they are fully due to scattering and obey the f^4 law.
          do i=1,size(signalOut%frequencies)
            do maf=1, quantity%template%noInstances
              do mif=1, quantity%template%noSurfs
                if ( iand ( ichar ( Quantity%mask( i + &
                  & (mif-1) * size(signalOut%frequencies), maf)), m_cloud) /= 0 ) then
                  quantity%values(i+(mif-1)*size(signalOut%frequencies), maf ) = &
                    & sourceQuantity%values ( myChannel + &
                    & (mif-1)*size(signalIn%frequencies), maf) *freq2(i)**4/ &
                    & (lsbFraction%values ( myChannel, 1 ) * freqL1(MyChannel)**4 + &
                    & usbFraction%values ( myChannel, 1 ) * freqU1(MyChannel)**4)
                end if
              end do
            end do
          end do
        else
          ! If both sidebands are within 20GHz but have very different penetration depths,
          ! where one is optically thick and one is optically thin. We may use a reference
          ! cloud radiance near the optically-thin sideband to split cloud radiance in the
          ! optically-thick sideband. The reference radiance must be a single sideband
          ! radiance in this case, namely usb (upper sideband in most cases). The usb is
          ! scaled to the sourceQty upperside (thin) via the f^4 law, and the rest cloud
          ! radiance is for the lower sideband.
          !
          ! In this case, channel will for the reference signal and sideband frac remains for
          ! source signal.
          !
          if ( .not. ValidateVectorQuantity ( usb, &
            & quantityType=(/l_cloudInducedRadiance/), sideband=(/-1,1/), &
            & minorFrame=.true. )) call Announce_Error ( key, no_error_code, &
            & 'Inappropriate reference radiance for splitting' )

          signalRef = GetSignal ( usb%template%signal )

          call allocate_test ( freq, size(signalRef%frequencies), 'frequencies', ModuleName )
          call allocate_test ( freqL, size(signalRef%frequencies), &
            & 'LSBfrequencies', ModuleName )
          call allocate_test ( freqU, size(signalRef%frequencies), &
            & 'USBfrequencies', ModuleName )
          freq = signalRef%centerFrequency + signalRef%direction*signalRef%frequencies
          freqL = signalOut%lo - freq    ! lower sideband freq
          freqU = signalOut%lo + freq    ! upper sideband freq
          if ( usb%template%sideband == -1) freq = freqL
          if ( usb%template%sideband == 1) freq = freqU
          do i=1,size(signalOut%frequencies)

            ! need to scale the opposite sideband for the output signal
            if ( quantity%template%sideband == 1 ) then
              ratio1 = lsbFraction%values(i,1)
              ratio2 = usbFraction%values(i,1)
              freq2 = freqL2
            end if
            if ( quantity%template%sideband == -1 ) then
              ratio2 = lsbFraction%values(i,1)
              ratio1 = usbFraction%values(i,1)
              freq2 = freqU2
            end if

            do maf=1, quantity%template%noInstances
              do mif=1, quantity%template%noSurfs
                scaledRad = usb%values(MyChannel+(mif-1)*size(signalIn%frequencies), maf) * &
                  & freq2(i)**4/freq(MyChannel)**4

                if ( iand ( ichar ( Quantity%mask( i + &
                  & (mif-1)*size(signalOut%frequencies), maf)), m_cloud) /= 0 ) then
                  quantity%values(i+(mif-1)*size(signalOut%frequencies), maf) = &
                    &  (sourceQuantity%values(i+(mif-1)*size(signalIn%frequencies), maf) - &
                    &  ratio1*scaledRad)/ratio2
                end if
              end do
            end do
          end do

          call deallocate_test ( freq, 'frequencies', ModuleName )
          call deallocate_test ( freqL,'LSBfrequencies', ModuleName )
          call deallocate_test ( freqU,'USBfrequencies', ModuleName )
        end if
      end if

      call deallocate_test ( freq1, 'frequencies', ModuleName )
      call deallocate_test ( freqL1,'LSBfrequencies', ModuleName )
      call deallocate_test ( freqU1,'USBfrequencies', ModuleName )
      call deallocate_test ( freq2, 'frequencies', ModuleName )
      call deallocate_test ( freqL2,'LSBfrequencies', ModuleName )
      call deallocate_test ( freqU2,'USBfrequencies', ModuleName )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromSplitSideband' )
    end subroutine FromSplitSideband

    ! ------------------------------------- GPHPrecision ----
    subroutine GPHPrecision ( key, quantity, &
      & tempPrecisionQuantity, refGPHPrecisionQuantity )
      ! Fill the GPH precision from the temperature and refGPH precision,
      ! ignoring the of diagonal elements (not available outside
      ! RetrievalModule anyway).
      integer, intent(in) :: key          ! For messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: TEMPPRECISIONQUANTITY
      type (VectorValue_T), intent(in) :: REFGPHPRECISIONQUANTITY

      ! Local variables

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.GPHPrecision', key )

      select case ( quantity%template%quantityType )
      case ( l_gph )
        if ( (tempPrecisionQuantity%template%noSurfs /= &
          &   quantity%template%noSurfs) .or. &
          &  (refGPHPrecisionQuantity%template%noInstances /= &
          &   quantity%template%noInstances) .or. &
          &  (tempPrecisionQuantity%template%noInstances /= &
          &   quantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_gph failed first test" )
        else if ( (any(quantity%template%surfs /= tempPrecisionQuantity%template%surfs)) .or. &
          & (any(quantity%template%phi /= tempPrecisionQuantity%template%phi)) .or. &
          & (any(quantity%template%phi /= refGPHPrecisionQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            &  "case l_gph failed second test" )
        else
          call GetGPHPrecision ( tempPrecisionQuantity, refGPHPrecisionQuantity, quantity%values )
        end if
      case default
        call Announce_error ( 0, no_error_code, 'GPH precision needed for result of GPHPrecision' )
      end select

    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end( 'FillUtils_1.GPHPrecision' )

    end subroutine GPHPrecision

      ! ------------------------------------- IWCFromExtinction ----
    subroutine IWCFromExtinction ( quantity, &
     & sourceQuantity, temperatureQuantity)
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (IWC) Quantity to fill
      type (VectorValue_T), intent(in) :: sourceQuantity ! cloud extinction
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)

      ! local variables
      real (r8), dimension(Temperaturequantity%template%noSurfs,quantity%template%noInstances) :: temp
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: zt, yt
      real (r8), dimension(sourceQuantity%template%noSurfs) :: ze, ye
      real (r8), dimension(quantity%template%noSurfs) :: z, Tz, Ez, iwc0
      real (r8), dimension(quantity%template%noInstances) :: x2, y2
      real (r8), dimension(Temperaturequantity%template%noInstances) :: x1, y1
      integer :: i

      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.IWCFromExtinction' )
      call MLSMessageCalls( 'push', constantName='IWCFromExtinction' )
      if ( .not. (quantity%template%coherent .and. sourceQuantity%template%coherent &
         .and. Temperaturequantity%template%coherent)) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "one of the quantities is not coherent")
      if ( sourceQuantity%template%noInstances /= Quantity%template%noInstances) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "IWC and Extinction profile numbers are not matched")

      if ( TemperatureQuantity%template%noInstances /= Quantity%template%noInstances ) then
         x1 = TemperatureQuantity%template%phi(1,:)
         x2 = quantity%template%phi(1,:)
         do i=1,Temperaturequantity%template%noSurfs
           y1 = TemperatureQuantity%values(i, :)
           call InterpolateValues( x1, y1, x2, y2, method='Linear' )
           temp(i,:) = y2
         end do
      else
         temp = TemperatureQuantity%values
      end if

      if ( quantity%template%verticalCoordinate == l_pressure ) then
        z = -log10 ( quantity%template%surfs(:,1) )
      else
        z = quantity%template%surfs(:,1)
      end if

      if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then
        ze = -log10 ( sourceQuantity%template%surfs(:,1) )
      else
        ze = sourceQuantity%template%surfs(:,1)
      end if

      if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then
        zt = -log10 ( Temperaturequantity%template%surfs(:,1) )
      else
        zt = Temperaturequantity%template%surfs(:,1)
      end if

      do i=1, quantity%template%noInstances
        yt = Temp(:, i)
        call InterpolateValues( zt, yt, z,  Tz, method='Linear' )
        ye = sourceQuantity%values(:, i)*1000._r8  ! converted to 1/km
        call InterpolateValues( ze, ye, z, Ez, method='Linear' )
        ! see ATBD for the conversion based on the size distribution from
        ! McFarquhar and Heymsfield [1996]
        !
        iwc0 = 10**(-2.77+0.01*(Tz-273.15))
        where(Ez < 0._r8) Ez = 0._r8
        Ez = Ez**(1./1.4)

        quantity%values(:,i) = iwc0*Ez
      end do
      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.IWCFromExtinction' )

    end subroutine IWCFromExtinction

    ! ------------------------------------------- LOSVelocity ---
    subroutine LOSVelocity ( key, qty, tngtECI, scECI, scVel)
      ! A special fill from geometry arguments
      use Geometry, only: OMEGA => W
      integer, intent(in) :: KEY
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: TNGTECI
      type (VectorValue_T), intent(in) :: SCECI
      type (VectorValue_T), intent(in) :: SCVEL


      ! Local variables
      integer :: MAF                      ! Loop counter
      integer :: MIF                      ! Loop counter
      integer :: noMAFs                   ! Number of major frames
      integer :: noMIFs                   ! Number of minor frames for this module
      integer :: x,y,z                    ! Indicies into the vectors

      real (r8), dimension(3) :: tngtVel   ! Due to rotation of earth
      real (r8), dimension(3) :: los       ! Normalised line of sight vector

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.LOSVelocity', key )
      ! First check that things are OK.
      if ( .not. ValidateVectorQuantity ( qty, &
        & quantityType=(/l_losVel/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_none/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Quantity to fill is not a valid LOS Velocity' )
      if ( .not. ValidateVectorQuantity ( tngtECI, &
        & quantityType=(/l_tngtECI/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_xyz/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Tangent ECI quantity is not of an appropriate form' )
      if ( .not. ValidateVectorQuantity ( scECI, &
        & quantityType=(/l_scECI/), &
        & minorFrame=.true., &
        & frequencyCoordinate=(/l_xyz/) ) ) call Announce_Error ( key, No_Error_Code, &
        & 'Spacecraft ECI quantity is not of an approriate form' )
      if ( qty%template%instrumentModule /= tngtECI%template%instrumentModule ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'LOS Velocity and Tangent ECI quantities are not for the same module' )
      if ( .not. IsModuleSpacecraft ( scECI%template%instrumentModule ) ) &
        & call Announce_Error ( key, No_Error_Code, &
        & 'Spacecraft ECI quantity is not for the spacecraft' )

      noMAFs = qty%template%noInstances
      noMIFs = qty%template%noSurfs

      do maf = 1, noMAFs
        do mif = 1, noMIFs

          ! First compute the tangent point velocity in ECI coordinates due
          ! to the rotation of the earth.  This no doubt makes approximations
          ! due to the slight non alignment between the earth's rotation axis and
          ! the ECI z axis, but I'm going to ignore this.

          ! Work out the indices in 3*mif,maf space
          x = 1 + 3*(mif-1)
          y = x+1
          z = x+2

          tngtVel= omega* (/ -tngtECI%values(y,maf), &
            &                 tngtECI%values(x,maf), 0.0_r8 /)

          ! Now compute the line of sight direction normal
          los = tngtECI%values(x:z,maf) - scECI%values(x:z,maf)
          los = los / sqrt(sum(los**2))

          ! Now compute the net velocity in this direction.  For the moment I'll
          ! assume +ve means the sc and tp are moving apart, and -ve that they're
          ! getting closer.

          qty%values(mif,maf) = dot_product(tngtVel, los) - &
            &                   dot_product(scVel%values(x:z,maf), los)

          ! Note that even though x,y,z have been used up to now for a GHz/THz
          ! minor frame quantity, they're OK with this sc one too.
        end do
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.LOSVelocity' )
    end subroutine LOSVelocity

    ! ------------------------------------- NoRadsPerMIF -----
    subroutine NoRadsPerMif ( key, quantity, measQty, asPercentage )
      use BitStuff, only: BITEQ
      ! Count number of valid (i.e., not masked) radiances
      ! optionally compute it as a percentage of largest number possible
      ! The largest number possible takes into account
      ! np = the number with negative precisions 
      !      masked with bits m_linalg + m_ignore
      ! ns = the number subsetted "do not use"
      !      masked only with bit m_linalg
      ! nv = the number valid
      !      not masked with bit m_linalg (either alone or with others)
      ! ni = total number of instances
      ! Then ni = possible + ns
      ! and pct = 100 * nv/possible
      integer, intent(in) :: KEY
      type(VectorValue_T), intent(inout) :: QUANTITY
      type(VectorValue_T), intent(in) :: MEASQTY
      logical, intent(in), optional   :: asPercentage ! as % of 
      ! Local variables
      integer  :: possible
      integer  :: MIF, MAF               ! Loop counters
      integer  :: I0, I1                 ! Indices
      logical  :: pct

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.NoRadsPerMif', key )
      pct = .false.
      if ( present(asPercentage) ) pct = asPercentage
      possible = measQty%template%noChans
      ! Do some fairly limited checking.
      if ( .not. ValidateVectorQuantity ( measQty, quantityType=(/l_radiance/), &
        & signal=(/quantity%template%signal/), sideband=(/quantity%template%sideband/) ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Quantity and measurement quantity disagree' )
      if ( associated ( measQty%mask ) ) then
        do maf = 1, quantity%template%noInstances
          do mif = 1, quantity%template%noSurfs
            i0 = 1 +  ( mif-1 ) * measQty%template%noChans
            i1 = i0 + measQty%template%noChans - 1
            quantity%values ( mif, maf ) = count ( &
              & iand ( ichar ( measQty%mask ( i0:i1, maf ) ), M_LinAlg ) == 0 )
            possible = measQty%template%noChans - &
              & count( biteq( ichar(measQty%mask( i0:i1, maf )), M_linAlg) )
            possible = max(possible, 1)
            if ( pct ) &
              & quantity%values( mif, maf ) = 100*quantity%values( mif, maf )/possible
          end do
        end do
      else
        quantity%values = measQty%template%noChans
        if ( pct ) quantity%values = 100*quantity%values/possible
      end if
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.NoRadsPerMif' )
    end subroutine NoRadsPerMIF

    ! ------------------------------------ PhiTanWithRefraction --
    subroutine PhiTanWithRefraction ( key, quantity, &
      & h2o, orbIncline, ptan, refGPH, temperature, ignoreTemplate )

      use Geometry, only: EARTHRADA, EARTHRADB, GEODTOGEOCLAT
      use Hydrostatic_M, only: HYDROSTATIC
      use Phi_Refractive_Correction_m, only: PHI_REFRACTIVE_CORRECTION_UP
      use Refraction_m, only: REFRACTIVE_INDEX

      integer, intent(in) :: KEY          ! Tree node, for error messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! PhiTan quantity to update
      type (VectorValue_T), intent(in) :: H2O         ! Water vapor
      type (VectorValue_T), intent(in) :: OrbIncline  ! Orbital inclination
      type (VectorValue_T), intent(in) :: PTan        ! Tangent pressure
      type (VectorValue_T), intent(in) :: RefGPH      ! Reference GPH
      type (VectorValue_T), intent(in) :: TEMPERATURE ! Temperature
      logical, intent(in)              :: IGNORETEMPLATE

      real(rp), dimension(quantity%template%noInstances) :: CP2, CSQ, REQ, SP2
      real(rp) :: PhiCorrs(temperature%template%noInstances,temperature%template%noSurfs)
      real(rp), dimension(temperature%template%noInstances) :: REQS
      real(rv), dimension(temperature%template%noSurfs) :: Heights, N, PhiCorr, PS
      integer :: I, J ! Subscripts, loop inductors

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.PhiTanWithRefraction', key )
      call MLSMessageCalls( 'push', constantName='PhiTanWithRefraction' )
      ! More sanity checks
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%instrumentModule /= ptan%template%instrumentModule ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'PHITan and PTan quantities are not for the same module' )
        if ( any(shape(quantity%values)/=shape(ptan%values)) ) &
          & call Announce_Error ( key, No_Error_Code, &
          & 'PHITan and PTan quantities are not the same size' )
        if ( .not. ValidateVectorQuantity ( temperature, &
          & quantityType=(/l_temperature/), coherent=.true., stacked=.true., &
          & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
          & call Announce_error ( key, no_error_code, 'Problem with temperature quantity for phiTan fill' )
        if ( .not. ValidateVectorQuantity ( h2o, &
          & quantityType=(/l_vmr/), molecule=(/l_h2o/), coherent=.true., stacked=.true., &
          & frequencyCoordinate=(/l_none/), verticalCoordinate=(/l_zeta/) ) ) &
          & call Announce_error ( key, no_error_code, 'Problem with h2o quantity for phiTan fill' )
        if ( .not. ValidateVectorQuantity ( refGPH, &
          & quantityType = (/l_refGPH/), coherent=.true., stacked=.true., &
          & verticalCoordinate=(/l_zeta/), frequencyCoordinate=(/l_none/), noSurfs=(/1/) ) ) &
          & call Announce_Error ( key, badrefGPHQuantity )
      endif
      i = 0
      if ( .not. DoHGridsMatch ( temperature, refGPH ) ) i = i + 1
      if ( .not. DoHGridsMatch ( temperature, h2o ) ) i = i + 2
      if ( quantity%template%noInstances /= orbincline%template%noInstances ) i = i + 4
      if ( .not. DoVGridsMatch ( temperature, h2o ) ) i = i + 8
      if ( i /= 0 ) call Announce_Error ( key, no_error_code, &
        & ' coordinates for temperature/refGPH/H2O disagree', &
        & extraInfo = (/i/) )

      ! compute equivalent earth radius at uncorrected tangent phi nearest surface
      csq = (earthrada * earthradb)**2 / &
            & ((earthrada**2-earthradb**2)*sin(orbincline%values(1,:)*deg2rad)**2 + &
            &   earthradb**2)
      cp2 = cos(quantity%values(1,:)*deg2rad)**2
      sp2 = 1.0_rp - cp2
      req = 0.001_rp*sqrt((earthrada**4*sp2 + csq**2*cp2) / &
                        & (earthrada**2*cp2 + csq*sp2))

      ! Interpolate REQ to temperature/h2o/refGPH hGrid
      call InterpolateValues ( quantity%values(1,:), req, &
        & temperature%template%phi(1,:), reqs, 'Linear', extrapolate='Constant' )

      ! Get pressures corresponding to Temperature etc
      ps = 10.0**(-temperature%template%surfs(:,1))

      ! OK, do the refraction calculation on Temperature's grids
      do i = 1, temperature%template%noInstances

        ! Get heights.  Temperature and RefGPH are on same hGrids.
        ! RefGPH is in meters, but Hydrostatic wants it in km.
        call Hydrostatic ( GeodToGeocLat ( temperature%template%geodLat(1,i) ), &
          & temperature%template%surfs(:,1), temperature%values(:,i), &
          & temperature%template%surfs(:,1), &
          & refGPH%template%surfs(1,1), 0.001*refGPH%values(1,i), &
          & heights )
        heights = heights + reqs(i)
        ! Get refractive indices.  Temperature and H2O are on same grids.
        call refractive_index ( ps, temperature%values(:,i), n, h2o%values(:,i) )
        do j = 1, temperature%template%noSurfs
          call phi_refractive_correction_up ( n(j:), heights(j:), phiCorr(j:) )
          phiCorrs(i,j) = phiCorr(temperature%template%noSurfs) * rad2deg
        end do
      end do

      ! Now interpolate Phi corrections onto the PhiTan grid and apply them.
      ! Its zeta grid is the value of PTan; its phi grid is its own value.

      do j = 1, size(quantity%values,2)
        call interpolateValues ( h2o%template%phi(1,:), quantity%values(:,j), &
          &                      h2o%template%surfs(:,1), ptan%values(:,j), &
          &                      phiCorrs, quantity%values(:,j), update=.true. )
      end do ! j

      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.PhiTanWithRefraction' )
    end subroutine PhiTanWithRefraction

      ! ------------------------------------- RHIFromOrToH2O ----
    subroutine RHIFromOrToH2O ( key, quantity, &
     & sourceQuantity, temperatureQuantity, &
     & dontMask, ignoreZero, ignoreNegative, interpolate, &
     & markUndefinedValues, invert )
      ! Convert h2o vmr to %RHI for all instances, channels, surfaces
      ! (See Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
      !  humidity measurement: Method and validation" Read et. al.
      !  J. Geoph. Res. Dec. 2001 (106) D23)

      !  Method:
      ! (1) straight convert--all quantities must have the same shape
      !     (strictly we assume they have _all_ the same geolocations)
      ! (2) interpolate--all quantities may have different shapes
      !     (the interpolation will be along the vertical coordinate only)
      !     I.e., for xQuantity (where x can be h2o or temperature)
      !     if NoChans(xQuantity) /= NoChans(Quantity)
      !        => use only xQuantity(channel==1)
      !     if NoInstances(xQuantity) /= NoInstances(Quantity)
      !        => use only xQuantity(instance==1)
      !
      ! (3) if invert is TRUE, like (1) but its inverse: %RHI to h2o vmr
      integer, intent(in) :: key          ! For messages
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (rhi) Quantity to fill
      type (VectorValue_T), intent(in) :: sourceQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)
      logical, intent(in)           ::    dontMask    ! Use even masked values
      logical, intent(in)           ::    ignoreZero  ! Ignore 0 values of h2o
      logical, intent(in)           ::    ignoreNegative  ! Ignore <0 values
      logical, intent(in)           ::    interpolate ! If VGrids or HGrids differ
      logical, intent(in)           ::    markUndefinedValues ! as UNDEFINED_VALUE
      logical, intent(in)           ::    invert      ! %RHI -> vmr if TRUE

      ! Local variables
      integer ::                          Channel     ! Channel loop counter
      integer ::                          Chan_h2o    ! Channel loop counter
      integer ::                          Chan_T      ! Channel loop counter
      logical, parameter ::               DEEBUG_RHI = .false.
      integer                          :: dim
      integer ::                          I           ! Instances
      integer ::                          I_H2O       ! Instance num for values
      integer ::                          I_T         ! Instance num for values
      integer ::                          invs        ! 1 if invert, else -1
      integer ::                          QINDEX
      integer ::                          N           ! Num. of summed values
      logical                          :: matched_h2o_channels
      logical                          :: matched_h2o_instances
      logical                          :: matched_sizes
      logical                          :: matched_surfs
      logical                          :: matched_T_channels
      logical                          :: matched_T_instances
      integer ::                          S           ! Surface loop counter
      integer ::                          S_H2O       ! Instance num for surfs
      integer ::                          S_RHI       ! Instance num for surfs
      integer ::                          S_T         ! Instance num for surfs
      logical ::                          skipMe
      real (r8) ::                        T
      character(len=*), parameter ::      VMR_UNITS = 'vmr'
      integer ::                          VMR_UNIT_CNV
      logical ::                          wereAnySkipped
      ! These automatic arrays could cause trouble later
      ! You may consider declaring them as pointers and
      ! calling allocate_test and deallocate_test
      real (r8), dimension(quantity%template%noSurfs) :: &
       &                                  zeta, TofZeta, H2OofZeta
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: &
       &                                  zetaTemperature, oldTemperature
      real (r8), dimension(sourceQuantity%template%noSurfs) :: &
       &                                  zetaH2o, oldH2o
      real (r8), dimension(quantity%template%noSurfs, quantity%template%noInstances) :: &
       &                                  values
      ! Executable statements
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.RHIFromOrToH2O', key )
      call MLSMessageCalls( 'push', constantName='RHIFromOrToH2O' )
      values = 0.
      ! Let any undefined values be so marked (but not necessarily masked)
      ! An exceptionally dubious step -- should remove this idea
      if ( markUndefinedValues ) Quantity%values = UNDEFINED_VALUE
      ! Will we convert %RHI to vmr?
      if ( invert ) then
        invs = 1
      else
        invs = -1
      end if
      ! Do we need to internally convert the vmr units?
      if ( VMR_UNITS == 'ppmv' ) then
        vmr_unit_cnv = 6
      else if ( VMR_UNITS == 'ppbv' ) then
        vmr_unit_cnv = 9
      else
        vmr_unit_cnv = 0
      end if
      ! Check that all is well
      if ( invert .and. interpolate ) then
       call Announce_Error ( key, No_Error_code, &
        & ' RHIFromOrToH2O unable to invert and interpolate simultaneously' )
       go to 9
      end if
      matched_sizes = .true.
      do dim=1, 2
        matched_sizes = matched_sizes .and. &
        & .not. any( size(Quantity%values,dim) /= &
        &(/ size(sourceQuantity%values,dim), &
        & size(temperatureQuantity%values,dim) /)&
        & )
      end do
      if ( .not. (matched_sizes .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Incompatible quantities in RHIFromOrToH2O--' //&
        & '(unless interpolating, all must have same shape)' )
       go to 9
      end if
      matched_surfs = .true.
      matched_surfs = matched_surfs .and. &
       & .not. any( Quantity%template%noSurfs /= &
       &(/ sourceQuantity%template%noSurfs, &
       & temperatureQuantity%template%noSurfs /)&
       & )
      if ( .not. (matched_surfs .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Different vertical coords in RHIFromOrToH2O--' //&
        & '(unless interpolating, all must be on the same VGrid)' )
       go to 9
      end if
      matched_h2o_channels = &
       &   (sourceQuantity%template%noChans == Quantity%template%noChans)
      matched_h2o_instances = &
       &   (sourceQuantity%template%noInstances == Quantity%template%noInstances)
      matched_T_channels = &
       &   (TemperatureQuantity%template%noChans == Quantity%template%noChans)
      matched_T_instances = &
       &   (TemperatureQuantity%template%noInstances == Quantity%template%noInstances)
      wereAnySkipped = .false.
      ! Now let's do the actual conversion
      do i=1, quantity%template%noInstances
        if ( quantity%template%coherent ) then
          s_rhi = 1
        else
          s_rhi = i
        end if
        if ( sourceQuantity%template%coherent ) then
          s_h2o = 1
        else
          s_h2o = i
        end if
        if ( temperaturequantity%template%coherent ) then
          s_t = 1
        else
          s_t = i
        end if
        ! zeta must be in log(hPa) units
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          zeta = -log10 ( quantity%template%surfs(:,s_rhi) )
        else
          zeta = quantity%template%surfs(:,s_rhi)
        end if
        if ( interpolate .and. .not. matched_h2o_instances ) then
          i_h2o = 1
        else
          i_h2o = i
        end if
        if ( interpolate .and. .not. matched_T_instances ) then
          i_T = 1
        else
          i_T = i
        end if
        if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then
          zetah2o = -log10 ( sourceQuantity%template%surfs(:,s_h2o) )
        else
          zetah2o = sourceQuantity%template%surfs(:,s_h2o)
        end if
        if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then
          zetaTemperature = -log10 ( Temperaturequantity%template%surfs(:,s_T) )
        else
          zetaTemperature = Temperaturequantity%template%surfs(:,s_T)
        end if
        N = 0
        do Channel=1, quantity%template%noChans
          if ( interpolate .and. .not. matched_h2o_channels ) then
            Chan_h2o = 1
          else
            Chan_h2o = Channel
          end if
          if ( interpolate .and. .not. matched_T_channels ) then
            Chan_T = 1
          else
            Chan_T = Channel
          end if
          if ( interpolate ) then
            do s=1, sourceQuantity%template%noSurfs
              qIndex = Chan_h2o + (s-1)*sourceQuantity%template%noChans
              oldH2o(s) = sourceQuantity%values(qIndex, i_h2o)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2o, oldH2o, &
             & zeta, H2OofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, Temperaturequantity%template%noSurfs
              qIndex = Chan_T + (s-1)*Temperaturequantity%template%noChans
              oldTemperature(s) = TemperatureQuantity%values(qIndex, i_T)
            end do
            call InterpolateValues( zetaTemperature, oldTemperature, &
             & zeta, TofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
          else
            do s=1, quantity%template%noSurfs
              N = N + 1
              qIndex = Channel + (s-1)*quantity%template%noChans
              H2OofZeta(s) = sourceQuantity%values(qIndex, i)
              TofZeta(s) = TemperatureQuantity%values(qIndex, i)
            end do
          end if
          do s=1, quantity%template%noSurfs
            N = N + 1
            qIndex = Channel + (s-1)*quantity%template%noChans
            skipMe = .false.
            skipMe = skipMe .or. &
            & .not. dontMask .and. ( &
            & (ignoreNegative .and. H2OofZeta(s) < 0.0 ) &
            & .or. (ignoreZero .and. H2OofZeta(s) == 0.0 ) &
            & )
            ! But skip no matter what else if temperature illegal
            skipMe = skipMe .or. TofZeta(s) <= 0.0
            if ( .not. skipMe ) then
              T = TofZeta(s)
              ! Quantity%values(qIndex, i) = &
              values(qIndex, i) = &
               & H2OofZeta(s) &
               & * &
               & RHIFromH2O_Factor(T, zeta(qIndex), vmr_unit_cnv, invert)
            end if
            wereAnySkipped = wereAnySkipped .or. skipMe
          end do
        end do
      end do
      if ( .not. associated ( quantity%mask ) ) then
        quantity%values = values
      else
        where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
          quantity%values = values
        end where
      end if
      if ( DEEBUG_RHI ) then
        call output('rhi Num. instances: ', advance='no')
        if ( invert ) then
          call output(sourceQuantity%template%noInstances, advance='yes')
        else
          call output(quantity%template%noInstances, advance='yes')
        end if
        call output('  size(surfs,2) ', advance='no')
        call output(size(quantity%template%surfs,2), advance='yes')
        call output('Were any rhi left undefined? ', advance='no')
        call output(wereAnySkipped, advance='yes')
        call dump(zeta, 'zeta(-log hPa)')
        do s=1, quantity%template%noSurfs
          if ( invert ) then
            zeta(s) = 1000000*Quantity%values(s,1)
          else
            zeta(s) = 1000000*sourceQuantity%values(s,1)
          end if
        end do
        call dump(zeta, 'H2O(ppmv)')
        call dump(TemperatureQuantity%values(:,1), 'Temperature(K)')
        if ( invert ) then
          call dump(sourceQuantity%values(:,1), 'RHI(%)')
        else
          call dump(Quantity%values(:,1), 'RHI(%)')
        end if
      end if
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.RHIFromOrToH2O' )
    end subroutine RHIFromOrToH2O
!MJF
    ! ------------------------------------- RHIPrecisionFromOrToH2O ----
    subroutine RHIPrecisionFromOrToH2O ( key, quantity, &
     & sourcePrecisionQuantity, tempPrecisionQuantity, sourceQuantity, &
     & temperatureQuantity, &
     & dontMask, ignoreZero, ignoreNegative, interpolate, &
     & markUndefinedValues, invert )
      ! For precisions:
      ! Convert h2o vmr to %RHI for all instances, channels, surfaces
      ! (See Eq. 9 from "UARS Microwave Limb Sounder upper tropospheric
      !  humidity measurement: Method and validation" Read et. al.
      !  J. Geoph. Res. Dec. 2001 (106) D23)

      !  Method:
      ! (1) straight convert--all quantities must have the same shape
      !     (strictly we assume they have _all_ the same geolocations)
      ! (2) interpolate--all quantities may have different shapes
      !     (the interpolation will be along the vertical coordinate only)
      !     I.e., for xQuantity (where x can be h2o or temperature)
      !     if NoChans(xQuantity) /= NoChans(Quantity)
      !        => use only xQuantity(channel==1)
      !     if NoInstances(xQuantity) /= NoInstances(Quantity)
      !        => use only xQuantity(instance==1)
      !
      ! (3) Don't know how to invert precisions
      !
      ! Note: If either of TempPrecision or H2OPrecision is negative, then
      ! set Precision to negative, too (if NEGATIVETOO is TRUE)
      integer, intent(in) :: key          ! For messages
      ! Actually, the meaning of the next two is reversed if invert is TRUE)
      type (VectorValue_T), intent(inout) :: QUANTITY ! (rhi) Quantity to fill
      type (VectorValue_T), intent(in) :: sourcePrecisionQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: tempPrecisionQuantity ! T(zeta)
      type (VectorValue_T), intent(in) :: sourceQuantity ! vmr (unless invert)
      type (VectorValue_T), intent(in) :: temperatureQuantity ! T(zeta)
      logical, intent(in)           ::    dontMask    ! Use even masked values
      logical, intent(in)           ::    ignoreZero  ! Ignore 0 values of h2o
      logical, intent(in)           ::    ignoreNegative  ! Ignore <0 values
      logical, intent(in)           ::    interpolate ! If VGrids or HGrids differ
      logical, intent(in)           ::    markUndefinedValues ! as UNDEFINED_VALUE
      logical, intent(in)           ::    invert      ! %RHI -> vmr if TRUE

      ! Local variables
      integer ::                          Channel     ! Channel loop counter
      integer ::                          Chan_h2oPrecision    ! Channel loop counter
      integer ::                          Chan_TPrecision      ! Channel loop counter
      integer ::                          Chan_h2o    ! Channel loop counter
      integer ::                          Chan_T      ! Channel loop counter
      logical, parameter ::               DEEBUG_RHI = .false.
      integer                          :: dim
      integer ::                          I           ! Instances
      integer ::                          I_H2OPrecision       ! Instance num for values
      integer ::                          I_TPrecision         ! Instance num for values
      integer ::                          I_H2O       ! Instance num for values
      integer ::                          I_T         ! Instance num for values
      logical                          :: matched_h2oPrecision_channels
      logical                          :: matched_h2oPrecision_instances
      logical                          :: matched_h2o_channels
      logical                          :: matched_h2o_instances
      logical                          :: matched_sizes
      logical                          :: matched_surfs
      logical                          :: matched_TPrecision_channels
      logical                          :: matched_TPrecision_instances
      logical                          :: matched_T_channels
      logical                          :: matched_T_instances
      integer ::                          N           ! Num. of summed values
      logical, parameter ::               NEGATIVETOO = .true.
      integer ::                          QINDEX
      real (r8) ::                        qty_precision
      integer ::                          S           ! Surface loop counter
      integer ::                          S_H2O       ! Instance num for surfs
      integer ::                          S_RHI       ! Instance num for surfs
      integer ::                          S_T         ! Instance num for surfs
      logical ::                          skipMe
      character(len=*), parameter ::      VMR_UNITS = 'vmr'
      integer ::                          VMR_UNIT_CNV
      logical ::                          wereAnySkipped
      ! These automatic arrays could cause trouble later
      ! You may consider declaring them as pointers and
      ! calling allocate_test and deallocate_test
      real (r8), dimension(quantity%template%noSurfs) :: &
       &                                  zeta, TPrecisionofZeta, sourcePrecisionofZeta, TofZeta, H2OofZeta
      real (r8), dimension(TempPrecisionquantity%template%noSurfs) :: &
       &                                  zetaTempPrecision, oldTempPrecision
      real (r8), dimension(sourcePrecisionQuantity%template%noSurfs) :: &
       &                                  zetaH2oPrecision, oldH2oPrecision
      real (r8), dimension(Temperaturequantity%template%noSurfs) :: &
       &                                  zetaTemperature, oldTemperature
      real (r8), dimension(sourceQuantity%template%noSurfs) :: &
       &                                  zetaH2o, oldH2o
      real (r8), dimension(quantity%template%noSurfs, quantity%template%noInstances) :: &
       &                                  values
      ! Executable statements
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.RHIPrecisionFromOrToH2O', key )
      call MLSMessageCalls( 'push', constantName='RHIPrecisionFromOrToH2O' )
      values = 0.
      ! Let any undefined values be so marked (but not necessarily masked)
      ! An exceptionally dubious step -- should remove this idea
      if ( markUndefinedValues ) Quantity%values = UNDEFINED_VALUE
      ! Will we convert %RHI to vmr?
      ! if ( invert ) then
      ! call Announce_Error ( key, No_Error_code, &
      !  & ' RHIPrecisionFromOrToH2O unable to invert' )
      ! go to 9
      ! end if
      ! Do we need to internally convert the vmr units?
      if ( VMR_UNITS == 'ppmv' ) then
        vmr_unit_cnv = 6
      else if ( VMR_UNITS == 'ppbv' ) then
        vmr_unit_cnv = 9
      else
        vmr_unit_cnv = 0
      end if
      ! Check that all is well
      if ( invert .and. interpolate ) then
       call Announce_Error ( key, No_Error_code, &
        & ' RHIPrecisionFromOrToH2O unable to invert and interpolate simultaneously' )
       go to 9
      end if
      matched_sizes = .true.
      do dim=1, 2
        matched_sizes = matched_sizes .and. &
        & .not. any( size(Quantity%values,dim) /= &
        &(/ size(sourcePrecisionQuantity%values,dim), &
        & size(tempPrecisionQuantity%values,dim), &
        & size(sourceQuantity%values,dim), &
        & size(temperatureQuantity%values,dim) /)&
        & )
      end do
      if ( .not. (matched_sizes .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Incompatible quantities in RHIPrecisionFromOrToH2O--' //&
        & '(unless interpolating, all must have same shape)' )
       go to 9
      end if
      matched_surfs = .true.
      matched_surfs = matched_surfs .and. &
       & .not. any( Quantity%template%noSurfs /= &
       &(/ sourcePrecisionQuantity%template%noSurfs, &
       & tempPrecisionQuantity%template%noSurfs, &
       & sourceQuantity%template%noSurfs, &
       & temperatureQuantity%template%noSurfs /)&
       & )
      if ( .not. (matched_surfs .or. interpolate) ) then
       call Announce_Error ( key, No_Error_code, &
        & 'Different vertical coords in RHIPrecisionFromOrToH2O--' //&
        & '(unless interpolating, all must be on the same VGrid)' )
       go to 9
      end if
      matched_h2oPrecision_channels = &
       &   (sourcePrecisionQuantity%template%noChans == Quantity%template%noChans)
      matched_h2oPrecision_instances = &
       &   (sourcePrecisionQuantity%template%noInstances == Quantity%template%noInstances)
      matched_TPrecision_channels = &
       &   (TempPrecisionQuantity%template%noChans == Quantity%template%noChans)
      matched_TPrecision_instances = &
       &   (TempPrecisionQuantity%template%noInstances == Quantity%template%noInstances)
      matched_h2o_channels = &
       &   (sourceQuantity%template%noChans == Quantity%template%noChans)
      matched_h2o_instances = &
       &   (sourceQuantity%template%noInstances == Quantity%template%noInstances)
      matched_T_channels = &
       &   (TemperatureQuantity%template%noChans == Quantity%template%noChans)
      matched_T_instances = &
       &   (TemperatureQuantity%template%noInstances == Quantity%template%noInstances)
      wereAnySkipped = .false.
      ! Now let's do the actual conversion
      do i=1, quantity%template%noInstances
        if ( quantity%template%coherent ) then
          s_rhi = 1
        else
          s_rhi = i
        end if
        if ( sourceQuantity%template%coherent ) then
          s_h2o = 1
        else
          s_h2o = i
        end if
        if ( temperaturequantity%template%coherent ) then
          s_t = 1
        else
          s_t = i
        end if
        ! zeta must be in log(hPa) units
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          zeta = -log10 ( quantity%template%surfs(:,s_rhi) )
        else
          zeta = quantity%template%surfs(:,s_rhi)
        end if
        if ( interpolate .and. .not. matched_h2oPrecision_instances ) then
          i_h2oPrecision = 1
        else
          i_h2oPrecision = i
        end if
        if ( interpolate .and. .not. matched_TPrecision_instances ) then
          i_TPrecision = 1
        else
          i_TPrecision = i
        end if
        if ( interpolate .and. .not. matched_h2o_instances ) then
          i_h2o = 1
        else
          i_h2o = i
        end if
        if ( interpolate .and. .not. matched_T_instances ) then
          i_T = 1
        else
          i_T = i
        end if
        if ( sourcePrecisionQuantity%template%verticalCoordinate == l_pressure ) then
          zetah2oPrecision = -log10 ( sourcePrecisionQuantity%template%surfs(:,s_h2o) )
        else
          zetah2oPrecision = sourcePrecisionQuantity%template%surfs(:,s_h2o)
        end if
        if ( TempPrecisionquantity%template%verticalCoordinate == l_pressure ) then
          zetaTempPrecision = -log10 ( TempPrecisionquantity%template%surfs(:,s_T) )
        else
          zetaTempPrecision = TempPrecisionquantity%template%surfs(:,s_T)
        end if
        if ( sourceQuantity%template%verticalCoordinate == l_pressure ) then
          zetah2o = -log10 ( sourceQuantity%template%surfs(:,s_h2o) )
        else
          zetah2o = sourceQuantity%template%surfs(:,s_h2o)
        end if
        if ( Temperaturequantity%template%verticalCoordinate == l_pressure ) then
          zetaTemperature = -log10 ( Temperaturequantity%template%surfs(:,s_T) )
        else
          zetaTemperature = Temperaturequantity%template%surfs(:,s_T)
        end if
        N = 0
        do Channel=1, quantity%template%noChans
          if ( interpolate .and. .not. matched_h2oPrecision_channels ) then
            Chan_h2oPrecision = 1
          else
            Chan_h2oPrecision = Channel
          end if
          if ( interpolate .and. .not. matched_TPrecision_channels ) then
            Chan_TPrecision = 1
          else
            Chan_TPrecision = Channel
          end if
          if ( interpolate .and. .not. matched_h2o_channels ) then
            Chan_h2o = 1
          else
            Chan_h2o = Channel
          end if
          if ( interpolate .and. .not. matched_T_channels ) then
            Chan_T = 1
          else
            Chan_T = Channel
          end if
          if ( interpolate ) then
            do s=1, sourcePrecisionQuantity%template%noSurfs
              qIndex = Chan_h2oPrecision + (s-1)*sourcePrecisionQuantity%template%noChans
              oldH2oPrecision(s) = sourcePrecisionQuantity%values(qIndex, i_h2oPrecision)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2oPrecision, oldH2oPrecision, &
             & zeta, sourcePrecisionofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, TempPrecisionquantity%template%noSurfs
              qIndex = Chan_TPrecision + (s-1)*TempPrecisionquantity%template%noChans
              oldTempPrecision(s) = TempPrecisionQuantity%values(qIndex, i_TPrecision)
            end do
            call InterpolateValues( zetaTempPrecision, oldTempPrecision, &
             & zeta, TPrecisionofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, sourceQuantity%template%noSurfs
              qIndex = Chan_h2o + (s-1)*sourceQuantity%template%noChans
              oldH2o(s) = sourceQuantity%values(qIndex, i_h2o)
            end do
            ! Know the following about the procedure we will call:
            ! First pair of args are old(X,Y), next pair are new(X,Y)
            ! We want newY(newX) via linear interp. w/o extrapolating
            ! and mark undefined values among oldY with UNDEFINED_VALUE
            call InterpolateValues( zetah2o, oldH2o, &
             & zeta, H2OofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
            do s=1, Temperaturequantity%template%noSurfs
              qIndex = Chan_T + (s-1)*Temperaturequantity%template%noChans
              oldTemperature(s) = TemperatureQuantity%values(qIndex, i_T)
            end do
            call InterpolateValues( zetaTemperature, oldTemperature, &
             & zeta, TofZeta, &
             & 'Linear', extrapolate='Constant', &
             & badValue=real(UNDEFINED_VALUE, r8), &
             & missingRegions=.TRUE. )
          else
            do s=1, quantity%template%noSurfs
              N = N + 1
              qIndex = Channel + (s-1)*quantity%template%noChans
              sourcePrecisionofZeta(s) = sourcePrecisionQuantity%values(qIndex, i)
              TPrecisionofZeta(s) = TempPrecisionQuantity%values(qIndex, i)
              H2OofZeta(s) = sourceQuantity%values(qIndex, i)
              TofZeta(s) = TemperatureQuantity%values(qIndex, i)
            end do
          end if
          do s=1, quantity%template%noSurfs
            N = N + 1
            qIndex = Channel + (s-1)*quantity%template%noChans
            skipMe = .false.
            skipMe = skipMe .or. &
            & .not. dontMask .and. ( &
            & (ignoreNegative .and. H2OofZeta(s) < 0.0 ) &
            & .or. (ignoreZero .and. H2OofZeta(s) == 0.0 ) &
            & )
            ! But skip no matter what else if temperature illegal
            skipMe = skipMe .or. TofZeta(s) <= 0.0
            if ( .not. skipMe ) then
              if ( invert ) then
                call H2OPrecFromRHi( H2OofZeta(s), &
                 & TofZeta(s), zeta(qIndex), vmr_unit_cnv, &
                 & sourcePrecisionofZeta(s), TPrecisionofZeta(s), &
                 & qty_precision, negativeToo )
              else
                call RHIPrecFromH2O( H2OofZeta(s), &
                 & TofZeta(s), zeta(qIndex), vmr_unit_cnv, &
                 & sourcePrecisionofZeta(s), TPrecisionofZeta(s), &
                 & qty_precision, negativeToo )
              endif
              ! Quantity%values(qIndex, i) = qty_precision
              values(qIndex, i) = qty_precision
            end if
            wereAnySkipped = wereAnySkipped .or. skipMe
          end do
        end do
      end do
      if ( .not. associated ( quantity%mask ) ) then
        quantity%values = values
      else
        where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
          quantity%values = values
        end where
      end if
      if ( DEEBUG_RHI ) then
        call output('rhi Num. instances: ', advance='no')
        if ( invert ) then
          call output(sourceQuantity%template%noInstances, advance='yes')
        else
          call output(quantity%template%noInstances, advance='yes')
        end if
        call output('  size(surfs,2) ', advance='no')
        call output(size(quantity%template%surfs,2), advance='yes')
        call output('Were any rhi left undefined? ', advance='no')
        call output(wereAnySkipped, advance='yes')
        call dump(zeta, 'zeta(-log hPa)')
        do s=1, quantity%template%noSurfs
          if ( invert ) then
            zeta(s) = 1000000*Quantity%values(s,1)
          else
            zeta(s) = 1000000*sourceQuantity%values(s,1)
          end if
        end do
        call dump(zeta, 'H2O(ppmv)')
        call dump(TemperatureQuantity%values(:,1), 'Temperature(K)')
        if ( invert ) then
          call dump(sourceQuantity%values(:,1), 'RHI(%)')
        else
          call dump(Quantity%values(:,1), 'RHI(%)')
        end if
      end if
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.RHIPrecisionFromOrToH2O' )
    end subroutine RHIPrecisionFromOrToH2O
!MJF
    ! ----------------------------------- FromASCIIFile --------
    subroutine FromAsciiFile ( key, quantity, filename, badRange )
      use IO_STUFF, only: GET_LUN
      use MACHINE, only: IO_ERROR
      use MOREMESSAGE, only: ONEMOREMESSAGE => MLSMESSAGE
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      integer, intent(in) :: FILENAME   ! ASCII filename to read from
      real(r8), dimension(2), optional, intent(in) :: BADRANGE ! Range for missing data
      ! Local variables
      integer :: LUN                    ! Unit number
      integer :: STATUS                 ! Flag from open/close/read etc.
      character(len=1024) :: FILENAMESTR
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromAsciiFile', key )
      call get_lun ( lun, msg=.false. )
      if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "No logical unit numbers available" )
      call Get_String ( filename, filenameStr, strip=.true. )
      open ( unit=lun, file=filenameStr, status='old', form='formatted', &
        & access='sequential', iostat=status )
      if ( status /= 0 ) then
        call io_Error ( "Unable to open ASCII input file ", status, filenameStr )
        call ONEMOREMESSAGE( MLSMSG_Error, ModuleName, 'Error opening ASCII file at %l', &
          & (/ source_ref(key) /) )
      end if
      read ( unit=lun, fmt=*, iostat=status ) quantity%values
      if ( status /= 0 ) then
        call io_Error ( "Unable to read ASCII input file ", status, filenameStr )
        call ONEMOREMESSAGE( MLSMSG_Error, ModuleName, 'Error reading ASCII file %l', &
          & (/ source_ref(key) /) )
      end if
      close ( unit=lun, iostat=status )
      if ( status /= 0 ) then
        call io_Error ( "Unable to close ASCII input file ", status, filenameStr )
        call ONEMOREMESSAGE( MLSMSG_Error, ModuleName, 'Error closing ASCII file at %l', &
          & (/ source_ref(key) /) )
      end if
      if ( present ( badRange ) ) then
        if ( .not. associated ( quantity%mask ) ) call CreateMask ( quantity )
        where ( quantity%values >= badRange(1) .and. &
          & quantity%values <= badRange(2) )
          quantity%mask = char(ior(ichar(quantity%mask),M_LinAlg))
        end where
      end if
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromAsciiFile' )
    end subroutine FromAsciiFile

    ! ------------------------------------------- FromInterpolatedQty
    subroutine FromInterpolatedQty ( qty, source, key, &
      & ignoreTemplate, ptan )
      use MLSNumerics, only: Interpolate_Regular_To_Irregular
      use VectorsModule, only: CreateVectorValue
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(in) :: SOURCE
      integer, intent(in) :: KEY
      logical, intent(in) :: IGNORETEMPLATE
      type (VectorValue_T), optional :: PTAN ! press. values

      ! Local variables
      real (r8), dimension(:), pointer :: oldSurfs, newSurfs
      real (r8), dimension(:,:), pointer :: newValues
      logical :: mySurfs, myNewValues
      integer :: instance
      integer :: status
      logical :: verbose

      ! Executable code
      verbose = ( switchDetail(switches, 'fill') > 0 ) ! -Sfill1
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromInterpolatedQty', key )
      call MLSMessageCalls( 'push', constantName='FromInterpolatedQty' )
      if ( .not. ignoreTemplate ) then
        if ( .not. DoQtysDescribeSameThing ( qty, source ) ) then
          call Announce_error ( key, no_error_code, &
            & 'Mismatch in quantities' )
          go to 9
        end if
        if ( .not. doHGridsMatch ( qty, source ) .and. .not. present(ptan) ) then
          call Announce_error ( key, no_error_code, &
            & 'Mismatch in horizontal grid' )
          go to 9
        end if
        if ( qty%template%noInstances /= source%template%noInstances ) then
          call Announce_error ( key, no_error_code, &
            & 'Mismatch in num of instances' )
          go to 9
        end if
        if ( .not. doFGridsMatch ( qty, source )  .and. .not. present(ptan) ) then
          call Announce_error ( key, no_error_code, &
            & 'Mismatch in frequency grid' )
          go to 9
        end if
      end if

      ! Two cases here, one where we have to interpolate vertically (has to be
      ! single channel quantity).  The other where we don't.  Most of the latter
      ! cases can be handled by the code that calls this routine.  The exception
      ! is when we've used the force option to e.g. copy one band into another.
      if ( .not. doVGridsMatch ( qty, source ) ) then
        ! This quantity needs vertical interpolation
        ! These checks are for cases the code can't (yet) handle,
        ! may add this functionality later.
        if ( present(ptan) ) then
          ! Very (too?) forgiving
        elseif ( qty%template%noChans /= 1) then
          call Announce_error ( key, no_error_code, &
            & 'Code cannot (yet?) interpolate multi channel quantities' )
          go to 9
        elseif ( .not. all ( (/ qty%template%coherent, source%template%coherent /) ) ) then
          call Announce_error ( key, no_error_code, &
            & 'Code cannot (yet?) interpolate incoherent quantities' )
          go to 9
        end if

        ! Work out vertical coordinate issues
        if ( present(ptan) ) then
          oldSurfs => source%template%surfs ( :, 1 )
          newSurfs => ptan%values ( :, 1 )
          mySurfs = .false.
        elseif ( qty%template%verticalCoordinate == l_pressure ) then
          nullify ( oldSurfs, newSurfs )
          call Allocate_test ( oldSurfs, source%template%noSurfs, 'oldSurfs', ModuleName )
          call Allocate_test ( newSurfs, qty%template%noSurfs, 'newSurfs', ModuleName )
          oldSurfs = -log10 ( source%template%surfs(:,1) )
          newSurfs = -log10 ( qty%template%surfs(:,1) )
          mySurfs = .true.
        else
          oldSurfs => source%template%surfs ( :, 1 )
          newSurfs => qty%template%surfs ( :, 1 )
          mySurfs = .false.
        end if

        ! Work out if we have to obey the mask
        myNewValues = .false.
        if ( associated ( qty%mask ) ) &
          & myNewValues = any ( iand ( ichar(qty%mask(:,:)), m_fill ) /= 0 )
        if ( myNewValues ) then
          nullify ( newValues )
          call Allocate_test ( newValues, qty%template%instanceLen, &
            & qty%template%noInstances, 'myNewValues', ModuleName )
        else
          newValues => qty%values
        end if

        ! OK, do the work
        if ( qty%template%logBasis ) then
          ! Whyd did you ever want to print this?
          if ( any( (/ source%values <= 0. /) ) .and. verbose ) then
            call MLSMessage( MLSMSG_Warning, ModuleName, &
              & 'log basis source values <= 0 in FromInterpolateQty', &
              & status=status )
            if ( status == 0 ) then
              call dump( qty%template )
              call dump( source%template )
            endif
          endif
          call InterpolateValues ( &
            & oldSurfs, log ( max ( source%values, sqrt(tiny(0.0_r8)) ) ), &
            & newSurfs, newValues, &
            & method='Linear', extrapolate='Constant' )
          newValues = exp ( newValues )
        else if( source%template%noInstances /= qty%template%noInstances ) then
          ! You have elected to interpolate when you have different numbers
          ! of horizontal instances.
          ! The only reasonable thing is to assume homogeneity of source.
          if ( present(ptan) ) then
            if ( source%template%coherent .and. qty%template%minorFrame ) then
              if ( source%template%noChans > 1 .or. &
                 & qty%template%noChans > 1 ) then
                call announce_error ( key, no_error_code, &
                  & 'Cannot interpolate regular-to-irregular with more than one channel' )
                go to 9
              end if
              if ( .not. associated(qty%values) ) &
                & call createVectorValue ( qty, 'qty', moduleName )
              call interpolate_regular_to_irregular ( &
                & source%template%phi(:,1), oldSurfs, source%values, &
                & qty%template%phi, ptan%values, qty%values )
            else
              ! If newSurfs is ptan, we can interpolate from the first column
              ! of source separately to each profile/MAF of qty.
              do instance=1, qty%template%noInstances
                call InterpolateValues ( &
                  & oldSurfs, source%values(:,1), &
                  & ptan%values(:,instance), newValues(:,instance), &
                  & method='Linear', extrapolate='Constant' )
              end do
            end if
          else
            ! Otherwise, all we can do is interpolate the first column and
            ! spread it.
            call InterpolateValues ( &
              & oldSurfs, source%values(:,1), &
              & newSurfs, newValues(:,1), &
              & method='Linear', extrapolate='Constant' )
            do instance=2, qty%template%noInstances
              newValues(:,instance) = newValues(:,1)
            end do
          end if
        else
          call InterpolateValues ( &
            & oldSurfs, source%values, &
            & newSurfs, newValues, &
            & method='Linear', extrapolate='Constant' )
        end if

        if ( myNewValues ) then
          where ( iand ( ichar(qty%mask(:,:)),m_fill) == 0 )
            qty%values = newValues
          end where
          call Deallocate_test ( newValues, 'myNewValues', ModuleName )
        end if

        ! Tidy up
        if ( mySurfs ) then
          call Deallocate_test ( oldSurfs, 'oldSurfs', ModuleName )
          call Deallocate_test ( newSurfs, 'newSurfs', ModuleName )
        end if
      else ! ------------------------
        ! No interpolation needed, more like the
        ! case handled in the calling code, except we're more lenient
        ! If we have a mask and we're going to obey it then do so

        ! There was a bug without consequences as originally coded:
        ! "qty" below was "quantity" while "source" was "sourceQuantity"
        ! As this subroutine was then an internal procedure to MLSL2Fill,
        ! and the formal args in the call were "Quantity" and "sourceQuantity"
        ! the correct items were actually being referred to
        if ( associated(qty%mask) ) then
          where ( iand ( ichar(qty%mask(:,:)), m_Fill ) == 0 )
            qty%values(:,:) = source%values(:,:)
          end where
        else ! Otherwise, just blindly copy
          qty%values = source%values
        end if
      end if
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromInterpolatedQty' )

    end subroutine FromInterpolatedQty

    !=============================== FromLosGrid ====
    subroutine FromLosGrid ( key, Qty, LOS, &
      & Ptan, Re, noFineGrid, extinction, errorCode )

      ! This is to fill a l2gp type of quantity with a los grid type of quantity.
      ! The los quantity is a vector quantity that has dimension of (s, mif, maf),
      ! where s is the path along los.

      ! Linear interpolation is used to fill l2gp grids and unfilled grids are
      ! marked with the baddata flag (-999.)

      ! Dummy arguments
      integer, intent(in) :: key          ! For messages
      type (VectorValue_T), intent(in) :: LOS ! Vector quantity to fill from
      type (VectorValue_T), intent(in) :: Ptan ! tangent pressure
      type (VectorValue_T), intent(in) :: Re ! Earth's radius
      type (VectorValue_T), intent(inout) :: QTY ! Quantity to fill
      integer, intent(in) :: noFineGrid   ! make finer sGrid with this number
      logical, intent(in) :: extinction  ! Flag for extinction fill and calculation
      integer, intent(out) :: errorCode ! Error code

      ! Local variables
      integer :: i, j, maf, mif                ! Loop counter
      integer :: maxZ, minZ                    ! pressure range indices of sGrid
      integer :: noMAFs,noMIFs,noDepths
      integer, dimension(qty%template%noSurfs,qty%template%noInstances) :: cnt
      real (r8), dimension(qty%template%noSurfs,qty%template%noInstances) :: out
      real (r8), dimension(qty%template%noSurfs) :: outZeta, phi_out, beta_out
      real (r8), dimension(los%template%noChans) :: x_in, y_in, sLevel
      real (r8), dimension(los%template%noSurfs) :: zt
      real (r8), dimension(los%template%noChans*noFineGrid) :: betaFine, TransFine, SFine
      real (r8), dimension(los%template%noChans, &
        & los%template%noSurfs,los%template%noInstances) :: beta
      real (r8) :: ds, ColTrans

      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromLosGrid', key )

      errorCode=0

      ! Make sure this quantity is appropriate
      !    if ( .not. ValidateVectorQuantity(qty, coherent=.TRUE., stacked=.TRUE., &
      !      & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) then
      !      call output ( " quantity vertical grid in FromLOSgrid is not valid")
      !    end if

      if ( qty%template%verticalCoordinate == l_pressure ) then
        outZeta = -log10 ( qty%template%surfs(:,1) )
      else
        outZeta = qty%template%surfs(:,1)
      end if

      noMAFs=los%template%noInstances
      noMIFs=los%template%noSurfs
      ! Now, we use frequency coordinate as sGrid along the path
      noDepths=los%template%noChans
      sLevel = los%template%frequencies

      ! the input losQty is the increment of cloud transmission function by default.
      ! it is converted to cloud extinction if extinction flag is on.
      if ( extinction ) then
        ! both sGrid and sFineGrid are expected to be evenly spaced at present
        ds = sLevel(2)-sLevel(1)
        do i=1,noDepths
          do j=1,noFineGrid
            Sfine(j+(i-1)*noFineGrid) = sLevel(i)+(j-1._r8)*ds/noFineGrid
          end do
        end do

        do maf=1,noMafs
          do mif=1,noMifs
            do i=1,noDepths
              ! convert increments to derivatives
              y_in(i) = los%values(i+(mif-1)*noDepths,maf)/ds
            end do
            call InterpolateValues(sLevel,y_in,sFine,TransFine,method='Linear')
            ! calculate column transmission function by integrating
            ! the derivatives on fine grid
            do i=1,noFineGrid*noDepths
              betaFine(i) = 0._r8
              colTrans=0._r8
              do j=1,i
                colTrans=colTrans + transFine(j)*ds/noFineGrid
              end do
              colTrans = 1._r8 - colTrans
              if ( colTrans > 0.02_r8) betaFine(i)= transFine(i)/colTrans
            end do
            ! interpolate betaFine back to the coarser sGrid
            call InterpolateValues(sFine,betaFine,sLevel,beta(:,mif,maf),method='Linear')

          end do
        end do

      else

        do maf=1,noMafs
          do mif=1,noMifs
            do i=1,noDepths
              beta(i,mif,maf)=los%values(i+(mif-1)*noDepths,maf)
            end do
          end do
        end do

      end if

      ! initialize quantity
      do j = 1, qty%template%noInstances
        do i = 1, qty%template%noSurfs
          qty%values(i,j)=qty%template%badValue
          cnt(i,j)=0
          out(i,j)=0._r8
        end do
      end do

      do maf=1,noMAFs
        zt = ptan%values(:,maf)   ! noChans=1 for ptan
        zt = (zt+3.)*16.                      ! converted to height in km
        do mif=1,noMIFs
          if ( ptan%values(mif,maf) .gt. -2.5) cycle ! for testing
          ! find altitude of each s grid
          x_in = sLevel**2/2./(re%values(1,maf)*0.001_r8 + zt(mif))
          ! converted to zeta
          x_in = x_in/16. + ptan%values(mif,maf)
          ! find minimum and maximum pressures indices in sGrid
          do i = 2,qty%template%noSurfs-1
            if ( ptan%values(mif,maf) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & ptan%values(mif,maf) > (outZeta(i)+outZeta(i-1))/2.) &
              & minZ = i
          end do
          if ( ptan%values(mif,maf) < (outZeta(1)+outZeta(2))/2.) minZ=1
          if ( ptan%values(mif,maf) > outZeta(qty%template%noSurfs)) cycle ! goto next mif

          do i = 2,qty%template%noSurfs-1
            if ( x_in(noDepths) < (outZeta(i)+outZeta(i+1))/2. .and. &
              & x_in(noDepths) > (outZeta(i)+outZeta(i-1))/2.) &
              & maxZ = i
          end do
          if ( x_in(noDepths) < (outZeta(1)+outZeta(2))/2.) cycle    ! goto next mif
          if ( x_in(noDepths) > outZeta(qty%template%noSurfs)) maxZ=qty%template%noSurfs

          ! get phi along path for each mif ( phi is in degree)
          y_in = los%template%phi(mif,maf) &
            & - atan(sLevel/(re%values(1,maf)*0.001_r8 + zt(mif)))*rad2deg
          ! interpolate phi onto standard vertical grids
          call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),phi_out(minZ:maxZ), &
            & method='Linear')
          ! interpolate quantity to standard vertical grids
          y_in = beta(:,mif,maf)
          call InterpolateValues(x_in,y_in,outZeta(minZ:maxZ),beta_out(minZ:maxZ), &
            & method='Linear')
          ! interpolate quantity to standard phi grids
          do i=minZ,maxZ
            do j = 2, qty%template%noInstances-1
              if ( phi_out(i) .lt. &
                & (qty%template%phi(1,j)+qty%template%phi(1,j+1))/2. &
                & .and. phi_out(i) .ge. &
                & (qty%template%phi(1,j-1)+qty%template%phi(1,j))/2. ) then
                out(i,j)=out(i,j) + beta_out(i)
                cnt(i,j)=cnt(i,j)+1       !  counter
              end if
            end do
          end do
        end do                            ! End surface loop
      end do                              ! End instance loop
      ! average all non-zero bins
      where (cnt > 0) qty%values = out/cnt
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromLosGrid' )

    end subroutine FromLosGrid

    ! --------------------------------------------- ByManipulation ---
    subroutine ByManipulation ( quantity, a, b, &
      & manipulation, key, ignoreTemplate, &
      & spreadflag, dimList, &
      & c )
      use MANIPULATIONUTILS, only: MANIPULATE
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), pointer :: A
      type (VectorValue_T), pointer :: B
      real(rv), optional            :: C  ! constant "c" in manipulation
      integer, intent(in) :: MANIPULATION
      integer, intent(in) :: KEY        ! Tree node
      logical, intent(in) :: ignoreTemplate ! If set throw caution to the wind
      ! The following args are important only for statistical functions
      logical, intent(in) :: SPREADFLAG ! If set spread across summed dimension
      character(len=*), intent(in)    :: DIMLIST ! E.g., 's' to shift surfaces, not chans
      ! Local parameters

      ! The 1 and 2-way manipulations must be entered exactly as shown
      ! Other more general manipulations are automatically recognized
      ! (with autoRecognizeGeneralExp)
      
      ! "More general" means free use of '+', '-', '*', '/' and appropriate
      ! nesting between '(' and ')' where appropriate.
      
      ! Why not do away with the 1-way and 2-way manipulations that can be
      ! done by the general manipulation? E.g., (a+b)/2 is easily 
      ! handled already, why make it a special 2-way?

      ! Isn't there a way to encapsulate the idiom
      ! if ( .not. associated ( quantity%mask ) ) then
      ! .   .   .
      !   end where
      ! endif
      ! So we can shorten this considerably?
      ! Probably need to create a temp array the same shape as quantity%values

      ! Not listed below but also available are the manipulations 
      ! 'a^c' and 'c^a' (also called 'a**c' and 'c**a')
      logical, parameter :: autoRecognizeGeneralExp = .true.
      integer, parameter :: MAXSTRLISTLENGTH = 128
      ! Local variables
      character (len=1) :: ABNAME
      type (VectorValue_T), pointer :: AORB
      real(rv) :: cc
      integer :: I
      logical :: MAPFUNCTION
      integer :: NUMWAYS ! 1 or 2
      character (len=128) :: MSTR
      logical :: OKSOFAR
      logical :: StatisticalFunction
      logical :: USESC
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ByManipulation', key )
      call MLSMessageCalls( 'push', constantName='ByManipulation' )
      ! Currently we have a rather brain dead approach to this, so
      ! check that what the user has asked for, we can supply.
      call get_string ( manipulation, mstr, strip=.true. )
      mstr = lowercase(mstr)
      
      usesC = present(c)
      ! StatisticalFunction = ( FindFirst( valid1WayManipulations, mstr ) > 7 )
      StatisticalFunction = any( &
        & indexes( &
        &   mstr, &
        &   (/ 'count ', 'min   ', 'max   ', 'mean  ', 'median', 'rms   ', 'stddev' /) &
        &   ) &
        &  > 1 )
      MapFunction = ( index(mstr, 'map') > 0 )
      
      numWays = 2 
      okSoFar = .true.
      do i = 1, numWays ! 2
        if ( i == 1 ) then
          if ( .not. associated(a) ) cycle
          aorb => a
          abName = 'a'
        else
          if ( .not. associated(b) ) cycle
          aorb => b
          abName = 'b'
        end if

        ! For minor frame quantities, check that we're on the same page
        if ( StatisticalFunction ) then
          ! We don't check for anything
        elseif ( MapFunction ) then
          ! We don't check for anything
        elseif ( .not. ignoreTemplate ) then
          if ( quantity%template%minorFrame ) then
            okSoFar = okSoFar .and. aorb%template%minorFrame .and. &
              & quantity%template%signal == aorb%template%signal .and. &
              & quantity%template%sideband == aorb%template%sideband .and. &
              & quantity%template%frequencyCoordinate == aorb%template%frequencyCoordinate
          elseif ( mstr == 'a+b' .or. mstr == 'a-b' ) then
            ! For a+/-b these quantities must share a template
            okSoFar = okSoFar .and. quantity%template%name == aorb%template%name
          end if
        else
          okSoFar = okSoFar &
            & .and. &
            & ( &
            & spreadflag .or. &
            & quantity%template%noInstances == aorb%template%noInstances &
            & .and. quantity%template%instanceLen == aorb%template%instanceLen &
            & )
        end if

        if ( .not. okSoFar ) then
          call Announce_Error ( key, no_error_code, &
            & abName // ' is not of the same (or close enough) type as quantity' )
          go to 9
        end if
      end do

      select case ( mstr )
      case ( 'a|b' )
        ! At some point we must treat this '|' operator more generally, like '+', ..
        if ( .not. associated ( quantity%mask ) ) then
          quantity%values = ior ( nint(a%values), nint(b%values) )
        else
          where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
            quantity%values = ior ( nint(a%values), nint(b%values) )
          end where
        end if
      case default
        ! This should be one of the cases which use the constant "c"
        if ( .not. present(c) ) then
          if ( autoRecognizeGeneralExp ) then
            cc = 0._rv
          else
            ! How did we get here?
            call Announce_Error ( key, no_error_code, &
              & 'trim(mstr) manipulation but no c supplied' )
            go to 9
          endif
        else
          cc = c
        endif
        call Manipulate( quantity, a, b, cc, mstr, &
          & spreadflag, dimList )
      end select
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ByManipulation' )

    end subroutine ByManipulation

    ! ----------------------------------------- FromL1B ----
    ! Fills a quantity that is stored as a dataset in either the
    ! l1boa file (lat, lons, etc.)
    ! l1brad file (radiances, radiance precisions, etc.)
    ! Optionally supply PrecisionQuantity when reading a radiance
    ! so that we mask radiances where the corresponding precisions
    ! are negative or masked themselves
    
    ! Naming conventions:
    ! l1boa datasets are named differently depending on hdfversion
    ! if 4, they are some way I have forgotten
    ! If 5, they are divided among 3 groups: GHz, THz, and spacecraft
    
    ! radiances are named the same independent of hdf version
    subroutine FromL1B ( root, quantity, chunk, filedatabase, &
      & isPrecision, suffix, PrecisionQuantity, BOMask )
      use BitStuff, only: NEGATIVEIFBITPATTERNSET
      integer, intent(in) :: root
      type (VectorValue_T), INTENT(INOUT) ::        QUANTITY
      type (MLSChunk_T), INTENT(IN) ::              CHUNK
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      logical, intent(in)               ::          ISPRECISION
      integer, intent(in)                        :: SUFFIX
      type (VectorValue_T), INTENT(IN), optional :: PRECISIONQUANTITY
      integer, intent(in), optional :: BOMask ! A pattern of bits--
                                              ! set prec. neg. if matched
      ! Local variables
      integer                               :: BO_error
      type (l1bData_T)                      :: BO_stat
      integer                               :: channel
      character (len=132)                   :: MODULENAMESTRING
      character (len=132)                   :: NAMESTRING
      integer                               :: FLAG, NOMAFS, maxMIFs
      type (l1bData_T)                      :: L1BDATA
      type (MLSFile_T), pointer             :: L1BFile
      type (MLSFile_T), pointer             :: L1BOAFile
      integer                               :: ROW, COLUMN
      integer                               :: myBOMask

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromL1B', root )
      myBOMask = 0
      if ( present(BOMask) ) myBOMask = BOMask
      ! print *, 'Filling vector quantity from l1b'
      ! L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      L1BOAFile => GetMLSFileByType(filedatabase, content='l1boa')
      ! fileID = L1BFile%FileID%f_id

      select case ( quantity%template%quantityType )
      case ( l_ECRtoFOV )
        call GetModuleName( quantity%template%instrumentModule, nameString )
        nameString = AssembleL1BQtyName('ECRtoFOV', L1BOAFile%HDFVersion, .TRUE., &
          & trim(nameString))
      case ( l_L1BMAFBaseline )
        call GetSignalName ( quantity%template%signal, nameString, &
          & sideband=quantity%template%sideband, noChannels=.TRUE. )
        nameString = AssembleL1BQtyName(nameString, L1BOAFile%HDFVersion, .FALSE.)
      case ( l_l1bMIF_TAI )
        if ( L1BOAFile%HDFVersion == HDFVERSION_5 ) then
          call GetModuleName ( quantity%template%instrumentModule, nameString )
          nameString = AssembleL1BQtyName ('MIF_TAI', L1BOAFile%HDFVersion, .FALSE., &
            & trim(nameString) )
        else ! ??? MIF_TAI is goofy in HDF4 files -- no sc, no tp, no GHz....
          nameString = 'MIF_TAI'
        end if
      case ( l_LosVel )
        call GetModuleName ( quantity%template%instrumentModule, nameString )
        nameString = AssembleL1BQtyName ('LosVel', L1BOAFile%HDFVersion, .TRUE., &
          & trim(nameString) )
      case ( l_orbitInclination )
        nameString = AssembleL1BQtyName('OrbIncl', L1BOAFile%HDFVersion, .FALSE., &
          & 'sc')
      case ( l_ptan )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('ptan', L1BOAFile%HDFVersion, .FALSE., &
          & trim(nameString))
      case ( l_radiance )
        call GetSignalName ( quantity%template%signal, nameString, &
          & sideband=quantity%template%sideband, noChannels=.TRUE. )
        nameString = AssembleL1BQtyName(nameString, L1BOAFile%HDFVersion, .FALSE.)
      case ( l_scECI )
        nameString = AssembleL1BQtyName('ECI', L1BOAFile%HDFVersion, .FALSE., 'sc')
      case ( l_scGeocAlt )
        nameString = AssembleL1BQtyName('GeocAlt', L1BOAFile%HDFVersion, .FALSE., &
          & 'sc')
      case ( l_scVel )
        nameString = AssembleL1BQtyName('Vel', L1BOAFile%HDFVersion, .FALSE., 'sc')
      case ( l_scVelECI )
        nameString = AssembleL1BQtyName('VelECI', L1BOAFile%HDFVersion, .FALSE., &
          & 'sc')
      case ( l_scVelECR )
        nameString = AssembleL1BQtyName('VelECR', L1BOAFile%HDFVersion, .FALSE., &
          & 'sc')
      case ( l_tngtECI )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('ECI', L1BOAFile%HDFVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeocAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('GeocAlt', L1BOAFile%HDFVersion, .TRUE., &
          & trim(nameString))
      case ( l_tngtGeodAlt )
        call GetModuleName( quantity%template%instrumentModule,nameString )
        nameString = AssembleL1BQtyName('GeodAlt', L1BOAFile%HDFVersion, .TRUE., &
          & trim(nameString))
      case default
        call Announce_Error ( root, cantFromL1B )
      end select

      ! Perhaps will need to read bright object status from l1bOA file
      if ( isPrecision .and. myBOMask /= 0 ) then
        call GetModuleName ( quantity%template%instrumentModule, moduleNameString )
        moduleNameString = AssembleL1BQtyName('BO_stat', L1BOAFile%HDFVersion, .TRUE., &
          & trim(moduleNameString))
        call ReadL1BData ( L1BOAFile, moduleNameString, BO_stat, noMAFs, &
          & flag=BO_error, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & NeverFail= .true., &
          & dontPad=DONTPAD )
      end if

      ! Possibly modify namestring based on whether a suffix field was supplied
      ! or a /isPrecision switch was set
      
      if ( suffix /= 0 ) then
        call Get_String ( suffix, &
          & nameString(len_trim(nameString)+1:), strip=.true. )
      elseif ( isPrecision ) then
        nameString = trim(nameString) // PRECISIONSUFFIX
      end if

      L1BFile => GetL1bFile(filedatabase, namestring)
      if (associated(L1BFile)) then
        call ReadL1BData ( L1BFile, nameString, l1bData, noMAFs, flag, &
          & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & NeverFail= .false., &
          & dontPad=DONTPAD )
      else
         flag = 1 ! Just a trick, so we reuse the printing of error message
      end if
      ! If it didn't exist in the not-a-radiance case, then we'll fail here.
      if ( flag /= 0 ) then
        if ( any(quantity%template%quantityType == &
          & (/ l_radiance, l_L1BMAFBaseline /) ) ) then
          ! This is the case where it's a radiance we're after and it's missing
          ! and we will allow this because we permit level 1 to omit certain
          ! bands when they are turned off
          quantity%values = DEFAULTUNDEFINEDVALUE ! -1.0
          do column=1, size(quantity%values(1,:))
            do row=1, size(quantity%values(:,1))
              call MaskVectorQty ( quantity, row, column, M_LinAlg )
            end do
          end do
        else
          call Announce_Error ( root, errorReadingL1B )
        endif
        go to 9
      end if
      if ( quantity%template%noInstances /= size ( l1bData%dpField, 3 ) .or. &
        &  quantity%template%instanceLen /= &
        &   size ( l1bData%dpField, 1 ) * size ( l1bData%dpField, 2 ) ) then
        call output ( 'Quantity shape:' )
        call output ( quantity%template%instanceLen )
        call output ( ' ( ' )
        call output ( quantity%template%noChans )
        call output ( ', ' )
        call output ( quantity%template%noSurfs )
        call output ( ' ), ' )
        call output ( quantity%template%noInstances, advance='yes' )
        call output ( 'L1B shape:' )
        call output ( size ( l1bData%dpField, 1 ) )
        call output ( ', ' )
        call output ( size ( l1bData%dpField, 2 ) )
        call output ( ', ' )
        call output ( size ( l1bData%dpField, 3 ), advance='yes' )
        call Announce_Error ( root, no_error_code, 'L1B data is wrong shape' )
        go to 9
      end if

      if ( isPrecision .and. myBOMask /= 0 .and. BO_error == 0 ) then
        noMAFs = size(l1bData%dpField, 3)
        maxMIFs = l1bData%maxMIFs
        if ( switchDetail(switches, 'glob') > 0 ) then ! e.g., 'glob1'
          call output ( 'Quantity shape:' )
          call output ( quantity%template%instanceLen )
          call output ( ' ( ' )
          call output ( quantity%template%noChans )
          call output ( ', ' )
          call output ( quantity%template%noSurfs )
          call output ( ' ), ' )
          call output ( quantity%template%noInstances, advance='yes' )
          call output ( 'L1B shape:' )
          call output ( size ( l1bData%dpField, 1 ) )
          call output ( ', ' )
          call output ( size ( l1bData%dpField, 2 ) )
          call output ( ', ' )
          call output ( size ( l1bData%dpField, 3 ), advance='yes' )
          call outputNamedValue( 'shape' // trim(namestring), shape(l1bData%dpField) )
          call outputNamedValue( 'shape(BO_stat)', shape(BO_stat%intField) )
          call outputNamedValue( 'noMAFs', noMAFs )
          call outputNamedValue( 'maxMIFs', maxMIFs )
          call outputNamedValue( 'noChans', quantity%template%noChans )
        endif
        if ( switchDetail(switches, 'glob') > 1 ) then ! e.g., 'glob2'
          call dump( l1bData%dpField(1,:,:), '(Before applying bright object mask)' )
        endif
        do channel = 1, quantity%template%noChans
        l1bData%dpField(channel,:,:) = &
          & NegativeIfBitPatternSet( l1bData%dpField(channel,:,:), &
          & BO_stat%intField(1, 1:maxMIFs, 1:noMAFs), myBOMask )
        enddo
        if ( switchDetail(switches, 'glob') > 1 ) &
          & call dump( l1bData%dpField(1,:,:), '(After applying bright object mask)' )
        call DeallocateL1BData(BO_stat)
      end if

      quantity%values = RESHAPE(l1bData%dpField, &
        & (/ quantity%template%instanceLen, quantity%template%noInstances /) )
      if ( isPrecision ) then
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            if ( quantity%values(row, column) < 0.d0 ) &
              & call MaskVectorQty(quantity, row, column, M_LinAlg)
          end do
        end do
      else if ( present(precisionQuantity) ) then
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            if ( isVectorQtyMasked(precisionQuantity, row, column, M_LinAlg) ) &
              & call MaskVectorQty(quantity, row, column, M_LinAlg)
          end do
        end do
        do column=1, size(quantity%values(1, :))
          do row=1, size(quantity%values(:, 1))
            if ( precisionQuantity%values(row, column) < 0.d0 ) then
              call MaskVectorQty(quantity, row, column, M_Ignore)
              call MaskVectorQty(quantity, row, column, M_LinAlg)
            endif
          end do
        end do
      end if

      if ( switchDetail(switches, 'l1b') > -1 ) call Dump( l1bData )
      call DeallocateL1BData(l1bData)
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromL1B' )
    end subroutine FromL1B

    ! ------------------------------------------- FromL2AUX --
    subroutine FromL2AUX ( qty, l2aux, errorCode )
      type ( VectorValue_T), intent(inout) :: QTY
      type ( L2AUXData_T), intent(in) :: L2AUX
      integer, intent(inout) :: ERRORCODE
      ! Local variables
      integer :: FIRSTPROFILE
      integer :: LASTPROFILE
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromL2AUX' )
      errorCode = 0
      ! Work out which profile in the l2aux this belongs to
      firstProfile = qty%template%instanceOffset - qty%template%noInstancesLowerOverlap
      lastProfile = firstProfile + qty%template%noInstances - 1
      ! In the case of minor/major frame quanties, while instanceOffset is
      ! Numbered from zero, as in L1B, our array starts from 1.
      if ( qty%template%minorFrame .or. qty%template%majorFrame ) then
        firstProfile = firstProfile + 1
        lastProfile = lastProfile + 1
      end if
      ! Check that the dimensions are appropriate
      if ( firstProfile < lbound ( l2aux%values, 3 ) ) then
        errorCode = CantFromL2AUX
        go to 9
      end if
      if ( lastProfile > ubound ( l2aux%values, 3 ) ) then
        errorCode = CantFromL2AUX
        go to 9
      end if
      if ( size ( l2aux%values, 1 ) /= qty%template%noChans .or. &
        &  size ( l2aux%values, 2 ) /= qty%template%noSurfs ) then
        errorCode = CantFromL2AUX
        go to 9
      end if
      ! Do the fill
      qty%values = reshape ( l2aux%values ( :, :,  &
        & firstProfile : lastProfile ), &
        & (/ qty%template%instanceLen, qty%template%noInstances /) )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromL2AUX' )
    end subroutine FromL2AUX

    ! --------------------------------------- UsingMagneticModel --
    subroutine UsingMagneticModel ( qty, gphQty, key )
      use Geometry, only: SECPERYEAR
      use IGRF_INT, only: FELDC, FELDCOF, TO_CART
      type (VectorValue_T), intent(inout) :: QTY
      type (VectorValue_T), intent(inout) :: GPHQTY
      integer, intent(in) :: KEY
      ! Local variables
      real :: B(3)                      ! Magnetic field
      integer :: INSTANCE               ! Loop counter
      character(len=8) :: options
      integer :: SURF                   ! Loop counter
      integer :: SURFOR1                ! Index
      real    :: XYZ(3)                 ! lat, lon, height for to_cart

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.UsingMagneticModel', key )
      options = ' '
      if ( switchDetail(switches,'clean') > -1 ) options = '-c'
      if ( .not. ValidateVectorQuantity ( qty, quantityType=(/l_magneticField/), &
        & frequencyCoordinate=(/ l_xyz /) ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Quantity does not describe magnetic field' )
        go to 9
      end if
      if ( .not. ValidateVectorQuantity ( gphQty, quantityType=(/l_gph/), &
        & frequencyCoordinate=(/ l_none /), verticalCoordinate=(/l_zeta/) ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'GPH quantity does not describe gph field' )
        go to 9
      end if
      if ( .not. DoHGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Quantity and GPHQuanity do not share the same horizontal basis' )
        go to 9
      end if
      if ( .not. DoVGridsMatch ( qty, gphQty ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Quantity and GPHQuantity do not share the same vertical basis' )
        go to 9
      end if

      ! Assume the time is close enough to constant that one call to
      ! FELDCOF is accurate enough.

      call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

      do instance = 1, qty%template%noInstances
        do surf = 1, qty%template%noSurfs
          if ( qty%template%stacked ) then
            surfOr1 = 1
          else
            surfOr1 = surf
          end if
          ! Convert (/ lat(deg), lon(deg), height(km) /) to cartesian (e-radius)
          call to_cart ( real( (/ qty%template%geodLat(surfOr1,instance), &
            &                     qty%template%lon(surfOr1,instance), &
            &                     gphQty%values(surf,instance)*1.0e-3 /) ), xyz )
          ! Compute the field at and w.r.t. cartesian coordinates
          call feldc ( xyz, b )
          qty%values ( surf*3-2 : surf*3, instance) = b
        end do
      end do

      if ( switchDetail(switches,'mag') > -1 ) &
        & call dump ( qty, options=options )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.UsingMagneticModel' )

    end subroutine UsingMagneticModel

    ! ------------------------------------- Hydrostatically ----
    subroutine Hydrostatically ( key, quantity, &
      & temperatureQuantity, refGPHQuantity, h2oQuantity, &
      & orbitInclinationQuantity, phiTanQuantity, geocAltitudeQuantity, maxIterations, &
      & phiWindow, phiWindowUnits, chunkNo )
    use MANIPULATEVECTORQUANTITIES, only: FINDCLOSESTINSTANCES
      ! Various hydrostatic fill operations
      integer, intent(in) :: key          ! For messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: TEMPERATUREQUANTITY
      type (VectorValue_T), intent(in) :: REFGPHQUANTITY
      type (VectorValue_T), pointer :: H2OQUANTITY
      type (VectorValue_T), pointer :: ORBITINCLINATIONQUANTITY
      type (VectorValue_T), pointer :: PHITANQUANTITY
      type (VectorValue_T), pointer :: GEOCALTITUDEQUANTITY
      integer, intent(in) :: MAXITERATIONS
      real(r8), intent(in) :: PHIWINDOW
      integer, intent(in) :: PHIWINDOWUNITS
      integer, intent(in), optional :: chunkNo
      ! H2OQuantity and GeocAltitudeQuantity have to be pointers
      ! as they may be absent.

      ! Local variables
      integer, dimension(:), pointer :: CLOSESTTEMPPROFILES
      logical :: verbose

      ! Executable code
      verbose = ( switchDetail(switches, 'fill') > -1 )

      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.Hydrostatically', key )

      select case ( quantity%template%quantityType )
      case ( l_gph )
        if ( (temperatureQuantity%template%noSurfs /= &
          &   quantity%template%noSurfs) .or. &
          &  (refGPHQuantity%template%noInstances /= &
          &   quantity%template%noInstances) .or. &
          &  (temperatureQuantity%template%noInstances /= &
          &   quantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_gph failed first test" )
          go to 9
        end if
        if ( (any(quantity%template%surfs /= temperatureQuantity%template%surfs)) .or. &
          & (any(quantity%template%phi /= temperatureQuantity%template%phi)) .or. &
          & (any(quantity%template%phi /= refGPHQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            &  "case l_gph failed second test" )
          go to 9
        end if
        call GetBasisGPH ( temperatureQuantity, refGPHQuantity, quantity%values )
      case ( l_ptan )
        if ( (temperatureQuantity%template%noInstances /= &
          &   refGPHquantity%template%noInstances) .or. &
          &  (temperatureQuantity%template%noInstances /= &
          &   h2oQuantity%template%noInstances) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed first test" )
          go to 9
        end if
        if ( (any(refGPHquantity%template%phi /= temperatureQuantity%template%phi)) .or. &
          & (any(h2oQuantity%template%phi /= temperatureQuantity%template%phi)) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed second test" )
          go to 9
        end if
        if ( (.not. ValidateVectorQuantity(quantity, minorFrame=.true.) ) .or. &
          &  (.not. ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ) .or. &
          &  (quantity%template%instrumentModule /= &
          &   geocAltitudeQuantity%template%instrumentModule) ) then
          call Announce_Error ( key, nonConformingHydrostatic, &
            & "case l_ptan failed third test" )
           call output('ValidateVectorQuantity(quantity, minorFrame=.true.) ')
           call blanks(3)
           call output( &
           & ValidateVectorQuantity(quantity, minorFrame=.true., sayWhyNot=.true.), &
           & advance='yes')

           call output('ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true.) ')
           call blanks(3)
           call output( &
           & ValidateVectorQuantity(geocAltitudeQuantity, minorFrame=.true., sayWhyNot=.true.), &
           & advance='yes')

           call output('quantity%template%instrumentModule ')
           call blanks(3)
           call output( &
           & quantity%template%instrumentModule, &
           & advance='yes')

           call output('geocAltitudeQuantity%template%instrumentModule ')
           call blanks(3)
           call output( &
           & geocAltitudeQuantity%template%instrumentModule, &
           & advance='yes')
           go to 9
        end if
        call Get2DHydrostaticTangentPressure ( quantity, temperatureQuantity,&
          & refGPHQuantity, h2oQuantity, orbitInclinationQuantity, &
          & phiTanQuantity, geocAltitudeQuantity, maxIterations, &
          & phiWindow, phiWindowUnits, chunkNo )
      case default
        call Announce_error ( 0, 0, 'No such fill yet' )
      end select

    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.Hydrostatically' )
      if ( .not. verbose ) return
      nullify( closestTempProfiles )
      call Allocate_Test ( closestTempProfiles, phitanquantity%template%noInstances, &
        "closestTempProfiles", ModuleName )
      call FindClosestInstances ( temperatureQuantity, phitanquantity, &
        & closestTempProfiles )
      call Deallocate_Test ( closestTempProfiles, "closestTempProfiles", ModuleName )
      
    end subroutine Hydrostatically

    ! -------------------------------------- FromIsotope -----------

    subroutine FromIsotope ( quantity, sourceQuantity, &
              & ratioQuantity )
      ! This routine fills one vector from another, given an appropriate
      ! isotope ratio.

      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY ! Quantity to take vmr from
      type (VectorValue_T), intent(in) :: RATIOQUANTITY ! Isotope ratio information

      ! Local variables
      real (r8) :: FACTOR                 ! Multiplier to apply to sourceQuantity

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromIsotope' )

      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate quantity for isotope fill")

      if ( .not. ValidateVectorQuantity ( sourceQuantity, &
        & quantityType=(/ l_vmr /), frequencyCoordinate=(/ l_none /) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Inappropriate source quantity for isotope fill")

      if ( .not. ValidateVectorQuantity ( ratioQuantity, &
        & quantityType=(/ l_isotopeRatio /), frequencyCoordinate=(/ l_none /), &
        & noInstances=(/1/), noSurfs=(/1/) ) ) &
        &   call MLSMessage ( MLSMSG_Error, ModuleName, &
        &     "Inappropriate form/quantity for isotope ratio")


      if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) .or. &
        &  .not. DoVGridsMatch ( quantity, sourceQuantity ) ) &
        &    call MLSMessage ( MLSMSG_Error, ModuleName, &
        &      "Quantity and source quantity don't match for isotope fill" )

      if ( quantity%template%molecule == sourceQuantity%template%molecule ) &
        & call MLSMessage( MLSMSG_Error, ModuleName, &
        &   "Source and quantity both describe same molecule in isotope fill")

      if ( ratioQuantity%template%molecule == quantity%template%molecule ) then
        ! Going from parent to isotope
        factor = ratioQuantity%values(1,1)
      else if ( ratioQuantity%template%molecule == sourceQuantity%template%molecule ) then
        ! Going from isotope to parent
        factor = 1.0/ratioQuantity%values(1,1)
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Unable to understand isotope fill" )
      end if

      quantity%values = sourceQuantity%values * factor
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromIsotope' )

    end subroutine FromIsotope

    ! ---------------------------------- WithEstdNoise ---
    subroutine WithEstNoise ( quantity, radiance, &
      & sysTemp, nbw, integrationTime )

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (VectorValue_T), intent(in) :: RADIANCE ! Radiances to use in calculation
      type (VectorValue_T), intent(in) :: SYSTEMP ! System temperature
      type (VectorValue_T), pointer :: NBW ! Noise bandwidth
      real(r8), intent(in) :: INTEGRATIONTIME ! Integration time in seconds

      ! Local variables
      integer :: C                        ! Channel loop counter
      integer :: S                        ! Surface loop counter
      integer :: I                        ! Index into first dimension of values

      real (r8), dimension(:), pointer :: WIDTH ! Channel widths in MHz

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithEstNoise' )

      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/), minorFrame=.true.) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid quantity for estimated noise fill.")

      if ( .not. ValidateVectorQuantity ( radiance, &
        & quantityType=(/l_radiance/), minorFrame=.true.) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid quantity for estimated noise fill.")

      if ( ( radiance%template%signal /= quantity%template%signal ) .or. &
        &  ( radiance%template%sideband /= quantity%template%sideband ) ) &
        & call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & "Quantity and radiances not same signal for estimated noise fill.")

      if ( ( systemp%template%signal /= quantity%template%signal ) .or. &
        &  ( systemp%template%sideband /= quantity%template%sideband ) ) &
        & call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & "Quantity and system temperature not same signal for estimated noise fill." )

      if ( .not. ValidateVectorQuantity ( &
        & sysTemp, &
        & quantityType=(/l_systemTemperature/), &
        & verticalCoordinate=(/l_none/), &
        & noInstances=(/1/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Invalid system temperature quantity for estimated noise fill.")

      if ( associated ( nbw ) ) then
        if ( .not. ValidateVectorQuantity ( &
          & nbw, &
          & quantityType=(/l_noiseBandwidth/), &
          & verticalCoordinate=(/l_none/), &
          & noInstances=(/1/) ) ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Invalid noise bandwidth quantity for estimated noise fill.")
        if ( ( nbw%template%signal /= quantity%template%signal ) .or. &
          &    ( nbw%template%sideband /= quantity%template%sideband ) ) &
          & call MLSMEssage ( MLSMSG_Error, ModuleName, &
          & "Quantity and noise bandwidth not same signal for estimated noise fill." )
        width => nbw%values(:,1)
      else
        width => signals(radiance%template%signal)%widths(:)
      end if

      i = 1
      do s = 1, quantity%template%noSurfs
        do c = 1, quantity%template%noChans
          quantity%values(i,:) = ( radiance%values(i,:) + sysTemp%values(c,1) ) / &
            & sqrt ( integrationTime * 1e6 * width(c) )
          i = i + 1
        end do
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithEstNoise' )

    end subroutine WithEstNoise

    ! ----------------------------------------- WithReflectorTemperature ---
    subroutine WithReflectorTemperature ( key, quantity, phiZero, termsNode )
      integer, intent(in) :: KEY         ! Tree node for messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to fill
      real(r8), intent(in) :: PHIZERO   ! Offset term
      integer, intent(in) :: TERMSNODE

      ! Local variables
      integer :: I                      ! Loop counter
      integer, DIMENSION(2) :: UNITASARRAY ! Unit for value given
      real (r8), DIMENSION(2) :: VALUEASARRAY ! Value give

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithReflectorTemperature', key )
      if ( quantity%template%quantityType /= l_reflTemp ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Inappropriate quantity for reflector temperature fill' )

      ! Loop over fourier terms, coefficients are son i+2
      do i = 0, nsons ( termsNode ) - 2
        call expr_check ( subtree(i+2,termsNode), unitAsArray, valueAsArray, &
          & (/PHYQ_Temperature/), unitsError )
        if ( unitsError ) call Announce_error ( termsNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Temperature/) )
        ! Add in this coefficient
        if ( i == 0 ) then
          quantity%values = valueAsArray(1)
        else if ( mod ( i, 2 ) == 1 ) then
          ! A sine term
          quantity%values(1,:) = quantity%values(1,:) + 2 * valueAsArray(1) * &
            & cos ( Deg2Rad * ((i+1)/2) * ( quantity%template%phi(1,:) - phiZero ) )
        else
          ! A cosine term
          quantity%values(1,:) = quantity%values(1,:) - 2 * valueAsArray(1) * &
            & sin ( Deg2Rad * ((i+1)/2) * ( quantity%template%phi(1,:) - phiZero ) )
        end if
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithReflectorTemperature' )

    end subroutine WithReflectorTemperature

    ! ----------------------------------------- WithReichlerWMOTP -------------
    subroutine WithReichlerWMOTP ( tpPres, temperature )
      
      use wmoTropopause, only: EXTRATROPICS, TWMO
      ! Implements the algorithm published in GRL
      ! Loosely called the "Reichler" algorithm
      ! Ideas the same as in WithWMOTropopause
      ! But implemented differently
      ! 
      type (VectorValue_T), intent(inout) :: TPPRES ! Result
      type (VectorValue_T), intent(in) :: TEMPERATURE
      ! Local variables
      integer :: instance
      integer :: invert
      real :: MISSINGVALUE
      integer :: nlev
      integer :: nvalid
      real, dimension( size(temperature%values, 1) ) :: p ! Pa
      real, parameter :: pliml = 65.*100 ! in Pa
      real, parameter :: plimlex = 65.*100 ! in Pa
      real, parameter :: plimu = 550.*100 ! in Pa
      real, dimension( size(temperature%values, 1) ) :: t
      real :: trp
      real, dimension(:), pointer :: xyTemp, xyPress
      logical :: alreadyDumped
      ! logical, parameter :: DEEBUG = .false.
      ! Executable
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithReichlerWMOTP' )
      nullify( xyTemp, xyPress )
      nlev = size(temperature%values, 1)
      MISSINGVALUE = REAL( DEFAULTUNDEFINEDVALUE )
      ! Loop over the instances
      tpPres%values = MISSINGVALUE
      if ( DEEBUG ) then
      call output( 'num levels', advance='no' )
      call output( nlev, advance='yes' )
      call output( 'num instances', advance='no' )
      call output( temperature%template%noInstances, advance='yes' )
      call output( 'shape of tpPress values', advance='no' )
      call output( shape(tpPres%values), advance='yes' )
      endif
      if ( -temperature%template%surfs(1,1) .gt.&
        & -temperature%template%surfs(2,1) ) then
        invert=1
        ! p = refGPH%values(nlev:1:-1,instance)*100.  ! hPa > Pa
        p = 10.0**(-temperature%template%surfs(nlev:1:-1,1))
        p = p * 100.  ! hPa > Pa
      else
        invert=0
        ! p=refGPH%values(:,instance)*100.         ! hPa > Pa
        p = 10.0**(-temperature%template%surfs(:,1)) 
        p = p * 100.  ! hPa > Pa
      end if
      instanceLoop: do instance = 1, temperature%template%noInstances
      alreadyDumped = .false.
        ! Check the temperatures are in a sensible range
        if ( any ( &
          & temperature%values(:,instance) < 10.0 .or. &
          & temperature%values(:,instance) > 1000.0 ) ) then
          call output('invalid temperatures in this instance', advance='yes' )
          cycle instanceLoop
        end if
        ! check vertical orientation of data
        ! (twmo expects ordered from top downward)
        ! but surfs are -1.0 * log10(press)
        if ( invert == 1 ) then
          t = temperature%values(nlev:1:-1,instance)
        else
          t = temperature%values(:,instance)
        end if
        where ( t < 0. .or. t > 100000. )
          t = MissingValue
        end where
        nvalid = count( .not. isFillValue(t) )
        if ( nvalid < 2 ) then
          if ( DEEBUG ) then
          call output('not enough valid temperatures in this instance', advance='yes' )
          call output(count( .not. isFillValue(t) ), advance='yes' )
          call output(count( isFillValue(t) ), advance='yes' )
          call dump(temperature%values(:,instance), 'temperature%values')
          endif
          cycle
        endif
        call Allocate_test (xyTemp, nvalid, 'xyTemp', ModuleName )
        call Allocate_test (xyPress, nvalid, 'xyPress', ModuleName )
        call RemoveFillValues( t, MISSINGVALUE, xyTemp, &
          & p, xyPress )
        call twmo(nvalid, xyTemp, xyPress, plimu, pliml, trp)
        if ( .not. alreadyDumped .and. DEEBUG ) then
           alreadyDumped = .true.
           call dump(p, 'p ')
           call dump(t, 't ')
           call dump(xyTemp, 'xyTemp')
           call dump(xyPress, 'xyPress')
           call output( 'plimu, pliml, trp', advance='no' )
           call output( (/ plimu, pliml, trp /), advance='yes' )
        endif
        ! Don't let tropopause sink too low in "extra tropics"
        if ( trp < plimlex .and. &
          & extraTropics(temperature%template%geodLat(1, instance) ) ) then
          if ( DEEBUG ) then
          call output( 'trp too low in extra tropics', advance='no' )
          call output( (/ plimlex, trp /), advance='yes' )
          endif
          trp = MISSINGVALUE
        endif
        if ( trp > 0. .and. trp < 100000000. ) tpPres%values(1, instance) = trp/100
        call Deallocate_test ( xyTemp, 'xyTemp', ModuleName )
        call Deallocate_test ( xyPress, 'xyPress', ModuleName )
      end do instanceLoop
      if ( DEEBUG ) call dump( tpPres%values, 'tpPres%values' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithReichlerWMOTP' )
    end subroutine WithReichlerWMOTP

    ! ----------------------------------------- WithWMOTropopause ------
    subroutine WithWMOTropopause ( tpPres, temperature, refGPH, grid )
      use Geometry, only: GEODTOGEOCLAT
      use Hydrostatic_M, only: HYDROSTATIC

      type (VectorValue_T), intent(inout) :: TPPRES ! Result
      type (VectorValue_T), intent(in) :: TEMPERATURE
      type (VectorValue_T), intent(in) :: REFGPH
      type (VGrid_T), intent(in) :: GRID

      ! Local variables
      real(rv), dimension(grid%noSurfs) :: TFINE ! Temperature on fine grid
      real(rv), dimension(grid%noSurfs) :: HFINE ! Temperature on fine grid
      real(rv), dimension(grid%noSurfs) :: DTDH ! d(tFine)/d(hFine)
      real(rv), dimension(grid%noSurfs) :: DUMMYT ! Extra arg. for interpolateValues
      integer :: I                      ! Instance counter
      integer :: S500                   ! 500mb surface?
      integer :: LOWCANDIDATE           ! Possbile tb but below 500hPa
      integer :: S                      ! Surface counter
      integer :: S0                     ! Surface to start looking from
      integer :: DS                     ! No. surfs less than 2km above s
      logical :: VALIDTP                ! Flag

      ! Now the rules for the WMO tropopause are:

      ! 1. The "first tropopause" is defined as the lowest height at which the
      ! lapse rate decreases to 2K per kilometer or less, provided also that
      ! the average lapse rate between this height and all higher altitudes
      ! within 2 kilometers does not exceed 2K per kilometer.

      ! 2. If, above the first tropopause, the average lapse rate between any
      ! height and all higher altitudes within a 1-kilometer intend exceeds 3K
      ! per kilometer, then another tropopause is defined by the same criteria
      ! as under 1 above. This second tropopause may be within or above the
      ! 1-kilometer layer.

      ! There are also two qualifying remarks attached to the selection
      ! criteria. They are as follows:

      ! 1. A height below the 500-mb level is not designated as a tropopause
      ! unless the sounding reaches the 200-mb level and the height is the only
      ! height satisfying the above definitions.

      ! 2. When the second or higher tropopauses are being determined, the
      ! 1-kilometer interval with an average lapse rate of 3K per kilometer can
      ! occur at any height above the conventional tropopause and not only at a
      ! height more than 2 kilometers above the first tropopause.

      ! Stricly melding the WMO tropopause definition with our own
      ! 'linear interpolation between tie points' mantra is problematic.
      ! If we are strict we'd only ever give the one of our basis points as
      ! the solution.  So we're going to interpolate to a higher resolution
      ! using a spline which pop's out the derivatives too.

      ! Note here we assume that the quantities supplied are all valid
      ! (ie describe the right quantities, and all on the same horizontal
      ! grid).

      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithWMOTropopause' )
      call MLSMessageCalls( 'push', constantName='WithWMOTropopause' )
      ! Loop over the instances
      tpPres%values = 0.0
      instanceLoop: do i = 1, temperature%template%noInstances
        ! Check the temperatures are in a sensible range
        if ( any ( &
          & temperature%values(:,i) < 10.0 .or. &
          & temperature%values(:,i) > 1000.0 ) ) then
          cycle instanceLoop
        end if
        ! Perhaps think later about keeping track of the 'top'.
        ! Currently we're going to assert that we have data extending to
        ! pressures lower than 200mb.

        ! Interpolate temperature onto 'finer' pressure grid
        call InterpolateValues ( &
          & temperature%template%surfs(:,1), temperature%values(:,i), &
          & grid%surfs(:,1), tFine, method='Spline' )
        ! Now get the height for this.
        call Hydrostatic ( GeodToGeocLat ( temperature%template%geodLat(1,i) ), &
          & grid%surfs(:,1), tFine, grid%surfs(:,1), &
          & refGPH%template%surfs(1,1), refGPH%values(1,i), hFine )
        ! Note while much of the software thinks in meters, the hydrostatic
        ! routine works in km.
        ! Now do another spline 'interpolation' to get dTdH
        call InterpolateValues ( hFine, tFine, hFine, dummyT, &
          & dYbyDx=dtdh, method='Spline' )
        ! Negate dTdH to turn it into lapse rate
        dTdH = - dTdH

        ! Locate the 500mb surface
        call Hunt ( grid%surfs(:,1), -log10(500.0_r8), s500 )

        ! Find the 'first' tropopause.  We're not going to bother
        ! with the 'second'
        lowCandidate = 0
        validTP = .false.
        s0 = 1
        tpHunt: do
          ! Find first place where lapse rate less than 2K/km
          s = FindFirst ( dtdh(s0:) < 2.0 ) + s0 - 1
          ! If not found such a place give up looking
          if ( s == s0 - 1 ) exit tpHunt
          ! Find last surface within 2km of this one
          ds = FindFirst ( hFine(s:) > hFine(s) + 2.0 ) - 1
          ! If not got data 2km above s0 (FindFirst gave 0, so ds=-1), give up
          ! Also, if data on such coarse resolution that the next point above
          ! s0 is more than 2km away (FindFirst gave 1 so ds=0), also give up
          if ( ds <= 0 ) exit tpHunt
          ! Must have mean lapse rate less than 2K/km within 2km of s.
          validTp = sum ( dTdH ( s : s + ds - 1 ) / ds ) < 2.0
          ! If this is below 500mb, keep an eye on it and keep looking
          if ( s < s500 .and. validTP ) then
            validTP = .false.
            lowCandidate = s
          end if
          if ( validTP ) exit tpHunt
          ! Otherwise keep looking higher up
          s0 = s + 1
          if ( s0 > grid%noSurfs ) exit tpHunt
        end do tpHunt

        ! Now pick up the pieces of that complex logic
        if ( .not. validTP .and. lowCandidate > 0 ) then
          s = lowCandidate
          validTP = .true.
        end if

        ! If we never found one, leave 0 in the result and move onto
        ! the next instance
        if ( .not. validTP ) cycle instanceLoop

        ! OK, our tropopause is below surface s
        ! Now we do some interpolation to get
        ! the value we really want.
        if ( s > 1 ) then
          tpPres%values(1,i) = 10.0 ** ( -( &
            & grid%surfs(s-1,1) + &
            & ( grid%surfs(s,1) - grid%surfs(s-1,1) ) * ( dTdH(s-1) - 2.0 ) / &
            & ( dTdH(s-1) - dTdH(s) ) ) )
        else
          tpPres%values(1,i) = 10.0 ** ( - grid%surfs(s,1) )
        end if
      end do instanceLoop
      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithWMOTropopause' )
    end subroutine WithWMOTropopause

    ! -------------------------------------------- QualityFromChisq --------
    ! Compute Quality as a function of chi^2
    ! namely, Q = scale / chi^2
    ! The computation is done at the first vertical surface
    ! unless heightNode is supplied
    subroutine QualityFromChisq ( key, quantity, sourceQuantity, &
      & scale, heightNode, ignoreTemplate )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in)    :: SOURCEQUANTITY ! Chisq like quantity on which it's based
      real(r8), intent(in)                 :: SCALE     ! A scale factor
      integer, intent(in)                  :: HEIGHTNODE ! What height to use
      logical, intent(in)                  :: IGNORETEMPLATE
      ! Local variables
      integer, dimension(2) :: UNITASARRAY ! From expr
      real(r8), dimension(2) :: VALUEASARRAY ! From expr
      real(r8) :: HEIGHT                ! The height to consider
      integer :: SURFACE                ! Surface index
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.QualityFromChisq', key )
      call MLSMessageCalls( 'push', constantName='QualityFromChisq' )
      ! Do some sanity checking
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%quantityType /= l_quality ) call Announce_error ( key, no_error_code, &
          & 'Quality quantity must be quality' )
        if ( sourceQuantity%template%quantityType /= l_chisqBinned ) call Announce_error ( &
          & key, no_error_code, 'sourceQuantity must be of type chisqBinned' )
        if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) ) call Announce_error ( &
          & key, no_error_code, 'quantity and sourceQuantity do not have matching hGrids' )
      endif
      ! Work out the height
      if ( heightNode /= 0 ) then
        if ( nsons ( heightNode ) /= 2 ) call Announce_Error ( key, no_error_code, &
          & 'Only one height can be supplied for quality fill' )
        call expr_check ( subtree(2,heightNode) , unitAsArray, valueAsArray, &
          & (/PHYQ_Pressure/), unitsError )
        if ( unitsError ) call Announce_error ( heightNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Pressure/) )
        height = - log10 ( valueAsArray(1) )
        call Hunt ( sourceQuantity%template%surfs(:,1), height, surface, nearest=.true. )
      else
        surface = 1
      end if
      quantity%values(1,:) = 0.0_r8
      where ( sourceQuantity%values(surface,:) /= 0.0_r8 )
        quantity%values(1,:) = scale / sourceQuantity%values(surface,:)
      end where
      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.QualityFromChisq' )
    end subroutine QualityFromChisq

    ! -------------------------------------------- StatusQuantity --------
    subroutine StatusQuantity ( key, quantity, sourceQuantity, statusValue, &
      & minValue, maxValue, heightNode, &
      & additional, exact, ignoreTemplate )
      integer, intent(in) :: KEY        ! Tree node
      type ( VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type ( VectorValue_T), intent(in) :: SOURCEQUANTITY ! Quantity on which it's based
      integer, intent(in) :: STATUSVALUE
      real(r8), intent(in) :: MINVALUE     ! A scale factor
      real(r8), intent(in) :: MAXVALUE     ! A scale factor
      integer, intent(in) :: HEIGHTNODE ! What height to compare at
      logical, intent(in) :: ADDITIONAL ! Is this an additional flag or a fresh start?
      logical, intent(in) :: EXACT ! Set status to exact statusValue , don't OR values
      logical, intent(in) :: IGNORETEMPLATE
      ! Local variables
      real(r8) :: HEIGHT                ! The height to consider
      integer :: SURFACE                ! Surface index
      integer, dimension(2) :: UNITASARRAY ! From expr
      real(r8), dimension(2) :: VALUEASARRAY ! From expr
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.StatusQuantity', key )
      call MLSMessageCalls( 'push', constantName='StatusQuantity' )
      ! Do some sanity checking
      if ( .not. ignoreTemplate ) then
        if ( quantity%template%quantityType /= l_status ) call Announce_error ( key, no_error_code, &
          & 'status quantity must be status' )
        if ( .not. DoHGridsMatch ( quantity, sourceQuantity ) ) &
          & call Announce_error ( &
          & key, no_error_code, 'quantity and sourceQuantity do not have matching hGrids' )
      endif
      ! Work out the height
      if ( heightNode /= 0 ) then
        if ( nsons ( heightNode ) /= 2 ) call Announce_Error ( key, no_error_code, &
          & 'Only one height can be supplied for status fill' )
        if ( sourceQuantity%template%verticalCoordinate /= l_zeta ) &
          & call Announce_Error ( key, no_error_code, 'Bad vertical coordinate for sourceQuantity' )
        call expr_check ( subtree(2,heightNode) , unitAsArray, valueAsArray, &
          & (/PHYQ_Pressure/), unitsError )
        if ( unitsError ) call Announce_error ( heightNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Pressure/) )
        height = - log10 ( valueAsArray(1) )
        call Hunt ( sourceQuantity%template%surfs(:,1), height, surface, nearest=.true. )
      else
        surface = 1
      end if
      if ( .not. additional ) quantity%values = 0.0_r8
      if ( size(sourceQuantity%values, 2) /= size(quantity%values, 2) ) then
        if ( sourceQuantity%values(surface,1) > maxValue .or. &
          & sourceQuantity%values(surface,1) < minValue ) then
          if ( exact ) then
            quantity%values(1,:) = statusValue
          else
            quantity%values(1,:) = ior ( nint ( quantity%values(1,1) ), statusValue )
          endif
        endif
      else
        ! quantity%values = iand ( nint ( quantity%values ), not ( statusValue ) )
        if ( exact ) then
          where ( sourceQuantity%values(surface,:) > maxValue .or. &
            & sourceQuantity%values(surface,:) < minValue )
              quantity%values(1,:) = statusValue
          end where
        else
          where ( sourceQuantity%values(surface,:) > maxValue .or. &
            & sourceQuantity%values(surface,:) < minValue )
              quantity%values(1,:) = ior ( nint ( quantity%values(1,:) ), statusValue )
          end where
        endif
      endif
      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.StatusQuantity' )
    end subroutine StatusQuantity

    ! ------------------------------------------ UsingLeastSquares -----
    subroutine UsingLeastSquares  ( key, Quantity, SourceQuantity, ptanQuantity, &
      & channel, method, scaleInstances, scaleRatio, scaleSurfs )
      ! This fills a coherent Quantity from a a typically incoherent
      ! SourceQuantity using a least-squares approximation to a first-order
      ! Taylor series.

      use HFTI_M, only: HFTI

      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), pointer :: PTANQUANTITY
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: METHOD
      real(r8), intent(inout) :: ScaleInstances, ScaleRatio, ScaleSurfs

      ! Local variables
      real(r8), dimension(:,:), pointer :: INSTS ! Instance coordinates for source
      real(r8), dimension(:,:), pointer :: LSMatrix ! Least Squares problem matrix
      real(r8), dimension(:), pointer :: RHS
      real(r8), dimension(:,:), pointer :: SOURCEHEIGHTS ! might be ptan.
      real(r8), dimension(:,:), pointer :: SURFS ! Surface coordinates for source
      real(r8) :: W, WMAX               ! Weight if method == l_lsWeighted
      integer :: QS, QI, SS, SI         ! Loop counters
      integer :: KRANK                  ! Rank of least-squares solution
      integer :: MYCHANNEL              ! Channel or 1
      integer :: NRows
      integer :: NSourceQuant  ! Number of values in SourceQuantity

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.UsingLeastSquares', key )

      ! Check the output quantity
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true. ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Quantity for least-squares fill must be stacked and coherent' )

      ! Also should have the condition:
      !  quantityType = (/ sourceQuantity%template%quantityType /)
      ! However, the code does not yet really support the ability to do that.

      ! Work out source vertical coordinate
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) then
        nullify ( sourceHeights )
        call Allocate_test ( sourceHeights, sourceQuantity%template%nosurfs, &
          & sourceQuantity%template%noinstances, 'sourceHeights', ModuleName )
        sourceHeights = ptanQuantity%values
        if ( quantity%template%verticalCoordinate /= l_zeta .and. WARNIFVERTCOORDNOTZETA ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinate in quantity to fill is not zeta' )
      else
        sourceHeights => sourceQuantity%template%surfs
        if ( sourceQuantity%template%verticalCoordinate /= quantity%template%verticalCoordinate ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinates in least-squares fill do not match' )
      end if

      ! Get the channel number
      if ( sourceQuantity%template%frequencyCoordinate /= l_none .and. &
        & channel == 0 ) then
        call Announce_Error ( key, no_error_code, &
          & 'Must supply channel for this least-squares fill' )
        go to 9
      end if
      myChannel = max(channel,1)

      ! Setup coordinate arrays
      nullify ( surfs, insts )
      call Allocate_test ( surfs, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'surfs', ModuleName )
      call Allocate_test ( insts, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'insts', ModuleName )

      ! Get the vertical coordinates, a function of instance for
      ! incoherent quantities
      do ss = 1, sourceQuantity%template%noSurfs
        if ( sourceQuantity%template%coherent ) then
          surfs(ss,:) = sourceHeights(ss,1)
        else
          surfs(ss,:) = sourceHeights(ss,:)
        end if
      end do

!      call dump ( surfs, 'surfs' )

      ! Get the horizontal coordinates, a function of height for unstacked quantities.
      do si = 1, sourceQuantity%template%noInstances
        if ( sourceQuantity%template%coherent ) then
          insts(:,si) = sourceQuantity%template%phi(1,si)
        else
          insts(:,si) = sourceQuantity%template%phi(:,si)
        end if
      end do

!      call dump ( insts, 'insts' )

      ! Allocate the arrays
      nullify ( lsMatrix, RHS )
      if ( sourceQuantity%template%regular ) then
        nSourceQuant = sourceQuantity%template%noInstances * &
          &            sourceQuantity%template%noSurfs
      else
        nSourceQuant = sourceQuantity%template%noInstances * &
          &                sourceQuantity%template%instanceLen
      end if
      call allocate_test ( lsMatrix, nSourceQuant, 4, 'LSMatrix', moduleName )
      call allocate_test ( RHS, nSourceQuant, 'RHS', moduleName )

      if ( method == l_lsWeighted ) then
        if ( scaleInstances < 0.0 ) scaleInstances = &
          ( quantity%template%phi(1,quantity%template%noInstances) - &
            quantity%template%phi(1,1) ) / &
          ( quantity%template%noInstances - 1 )
        if ( scaleSurfs < 0.0 ) scaleSurfs = &
          ( quantity%template%surfs(quantity%template%noSurfs,1) - &
            quantity%template%surfs(1,1) ) / &
          ( quantity%template%noSurfs - 1 )
        scaleSurfs = scaleSurfs * scaleRatio
      end if
      ! Now loop over the output quantity points and set up the least-squares
      ! problem.
      do qi = 1, quantity%template%noInstances
        do qs = 1, quantity%template%noSurfs
          nRows = 0
          wMax = -1.0
          do si = 1, sourceQuantity%template%noInstances
            do ss = 1, sourceQuantity%template%noSurfs
              if ( method == l_lslocal ) then
                if ( surfs(ss,si) < quantity%template%surfs(max(qs-1,1),1) ) cycle
                if ( surfs(ss,si) > quantity%template%surfs(min(qs+1,quantity%template%noSurfs),1) ) cycle
                if ( insts(ss,si) < quantity%template%phi(1,max(qi-1,1)) ) cycle
                if ( insts(ss,si) > quantity%template%phi(1,min(qi+1,quantity%template%noInstances)) ) cycle
              end if
              nRows = nRows + 1
              lsMatrix(nRows,:3) = (/ 1.0_r8, &
                & insts(ss,si) - quantity%template%phi(1,qi), &
                & surfs(ss,si) - quantity%template%surfs(qs,1) /)
              lsMatrix(nRows,4) = 0.5 * lsMatrix(nRows,2) * lsMatrix(nRows,3)
              rhs(nRows) = sourceQuantity%values(myChannel+(ss-1)*sourceQuantity%template%noChans,si)
              if ( method == l_lsWeighted ) then
                w = exp( - (lsMatrix(nRows,2)/scaleInstances)**2 &
                       & - (lsMatrix(nRows,3)/scaleSurfs)**2 )
                if ( w > wMax ) then
                  wMax = w
                else if ( w < sqrt(epsilon(w))*wMax ) then
                  nRows = nRows - 1
                  cycle
                end if
                lsMatrix(nRows,:) = lsMatrix(nRows,:) * w
                rhs(nRows) = rhs(nRows) * w
              end if
            end do
          end do
          if ( nRows < 4 ) then
            quantity%values(qs,qi) = quantity%template%badValue
          else
            ! Solve the least squares problem
            call hfti ( lsMatrix(:nRows,:), rhs(:nRows), krank=krank )
            if ( krank < 4 ) then
              quantity%values(qs,qi) = quantity%template%badValue
            else
              quantity%values(qs,qi) = rhs(1)
            end if
          end if
        end do
      end do

      ! Now tidy up
      call deallocate_test ( surfs,    'surfs',    moduleName )
      call deallocate_test ( insts,    'insts',    moduleName )
      call deallocate_test ( lsMatrix, 'LSMatrix', moduleName )
      call deallocate_test ( RHS,      'RHS',      moduleName )
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) &
        & call Deallocate_test ( sourceHeights, 'sourceHeights', ModuleName )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.UsingLeastSquares' )

    end subroutine UsingLeastSquares

    !=============================== FromGrid ============
    subroutine FromGrid(quantity, grid, allowMissing, errorCode)
      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (GriddedData_T), intent(inout) :: GRID ! Grid to fill it from
      ! Needs to be inout because we wrap it.
      logical, intent(in) :: ALLOWMISSING ! If set missing data in grid ok
      integer, intent(out) :: ERRORCODE   ! Error code (one of constants defined above)

      ! Local variables
      logical :: check
      logical :: DEEBUG
      integer :: instance,surf            ! Loop counter
      integer :: instIndex,surfIndex      ! Indices
      real(rv) :: newValue
      logical :: noGrid
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromGrid' )
      errorCode = 0
      DEEBUG = .false.
      ! DEEBUG = ( grid%quantityName == 'TEMPERATURE' .or. &
      !   & grid%description == 'Temperature' )
      check = (switchDetail(switches, 'cgrid') > -1)
      if ( check ) call outputNamedValue( 'check', check )
      noGrid = .not. associated(grid%field)
      if ( .not. noGrid ) noGrid = grid%empty
      if ( noGrid ) then
        ! Must allow this as missing gmao files are a possibility
        ! to be handled with grace and aplomb
        call MLSMessage ( MLSMSG_Warning, moduleName, &
          & 'No tropopause or whatever values in grid--filling with missing values' )
        quantity%values = REAL( DEFAULTUNDEFINEDVALUE ) ! grid%missingValue
        go to 9
      end if

      if ( quantity%template%verticalCoordinate /= l_zeta .and. &
        & quantity%template%noInstances > 1 .and. grid%noHeights > 1 ) then
        errorCode=NotZetaForGrid
        go to 9
      end if

      instIndex=1
      surfIndex=1

      ! Wrap the grid to be sure that we can interpolate it in longitude
      ! This will skip out if it's already been done.
      call WrapGriddedData ( grid )

      if ( DEEBUG ) then
        call dump( grid )
        call dump( quantity%template%surfs, 'zeta surfs' )
        call outputNamedValue( '(noSurfs,noInstances)', &
          & (/ quantity%template%noSurfs,quantity%template%noInstances /) )
      endif
      do instance = 1, quantity%template%noInstances
        if ( .not. quantity%template%stacked) instIndex=instance

        do surf = 1, quantity%template%noSurfs
          if ( .not. quantity%template%coherent) surfIndex=surf
          if ( DEEBUG ) then
            call outputNamedValue( 'interp pressure', 10.0**(-quantity%template%surfs(surf,instIndex)) )
            if ( 10.0**(-quantity%template%surfs(surf,instIndex)) < 1.E-6 ) then
              call outputNamedValue( '(surf,instIndex)', (/ surf,instIndex /) )
            endif
          endif
          call l3ascii_interp_field(grid, newValue, &
            & pressure=10.0**(-quantity%template%surfs(surf,instIndex)), &
            & lat=quantity%template%geodLat(surfIndex,instance), &
            & lon=quantity%template%lon(surfIndex,instance), &
            & lst=quantity%template%solarTime(surfIndex,instance), &
            & sza=quantity%template%solarZenith(surfIndex,instance), &
            & date=quantity%template%time(surfIndex,instance), &
            & debug=.false. )
            ! & debug=(instance==1 .and. surf==1 .and. grid%quantityName == 'CO') )
          if ( newValue >= nearest ( grid%missingValue, -1.0 ) .and. &
            &  newValue <= nearest ( grid%missingValue,  1.0 ) .and. &
            & .not. allowMissing ) errorCode = MissingDataInGrid
          quantity%values(surf,instance) = newValue
        end do                            ! End surface loop
      end do                              ! End instance loop
      
      if ( check ) then
          call l3ascii_interp_field(grid, newValue, &
            & pressure=100.0_rv, &
            & lat=0.0_rv, &
            & lon=0.0_rv, &
            & lst=0.0_rv, &
            & sza=0.0_rv, &
            & date=quantity%template%time(1,1), &
            & debug=.false. )
          call outputNamedValue( 'interpolated value', newValue )
      endif
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromGrid' )
    end subroutine FromGrid

    !=============================== FromL2GP ==========
    subroutine FromL2GP ( quantity,l2gp, interpolate, profile, &
      & errorCode, ignoreGeolocation, fromPrecision )

      ! If the times, pressures, and geolocations match, fill the quantity with
      ! the appropriate subset of profiles from the l2gp

      ! Dummy arguments
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      type (L2GPData_T), intent(in), target :: L2GP ! L2GP to fill from
      logical, intent(in) :: interpolate  ! Flag
      integer, intent(in) :: profile    ! Single profile to use or -1 for default
      integer, intent(out) :: errorCode ! Error code
      logical, intent(in) :: ignoreGeolocation  ! Flag
      logical, intent(in) :: fromPrecision ! Flag

      ! Local parameters
      real(r8), parameter :: FTOL = 1.0e-3 ! 1 kHz
      real(r8), parameter :: TOLERANCE=0.05 ! Tolerence for angles
      real(r8), parameter :: TIMETOL=5.0  ! Tolerence for time (not sure why this
      ! needs to be so big !????????? NJL)

      ! Local variables
      integer ::    FIRSTPROFILE, LASTPROFILE
      integer, dimension(1) :: FIRSTPROFILEASARRAY
      integer :: INSTANCE               ! Loop counter
      integer :: THISPROFILE            ! Index
      type (Coefficients_R8) :: COEFFS  ! For interpolation
      real (r8), dimension(quantity%template%noSurfs) :: outZeta
      real (r4), dimension(:,:,:), pointer :: SOURCE

      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromL2GP' )
      call MLSMessageCalls( 'push', constantName='FromL2GP' )
      errorCode=0
      ! Make sure this quantity is appropriate
      if ( .not. ValidateVectorQuantity(quantity, coherent=.TRUE., stacked=.TRUE., &
        & verticalCoordinate= (/ l_pressure, l_zeta /) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity to fill is not on pressure or zeta coordinates' )

      if ( (quantity%template%noChans/=l2gp%nFreqs) .and. &
        &  ((quantity%template%noChans/=1) .or. (l2gp%nFreqs/=0)) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of channels' )
!       call dump ( l2gp%frequency, 'l2gp' )
!       call dump ( quantity%template%frequencies, 'quantity' )
!       if ( associated ( quantity%template%frequencies ) ) then
!         if ( any ( abs ( l2gp%frequency - &
!           & quantity%template%frequencies ) > fTol ) ) &
!           & call MLSMessage ( MLSMSG_Error, ModuleName, &
!           & 'Quantity and l2gp have different frequency grids' )
!       end if

      if ( quantity%template%noSurfs /= l2gp%nLevels .and. (.not. interpolate) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and l2gp have different number of surfaces (set interpolate?)' )

      if ( .not. interpolate ) then
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          if ( any(ABS(-LOG10(quantity%template%surfs(:,1))+ &
            & LOG10(l2gp%pressures)) > TOLERANCE) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity and l2gp are on different surfaces (set interpolate?)' )
        else                                ! Must be l_zeta
          if ( any(ABS(quantity%template%surfs(:,1)+ &
            & LOG10(l2gp%pressures)) > TOLERANCE) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Quantity and l2gp are on different surfaces (set interpolate?)' )
        end if
      end if

      ! Skip the position checks if we're forcing in a particular profile.
      firstProfile = 1
      lastProfile=firstProfile+quantity%template%noInstances-1
      if ( (profile == -1) .and. (.not. ignoreGeolocation) ) then
        ! Attempt to match up the first location
        firstProfileAsArray=minloc(abs(quantity%template%phi(1,1)-l2gp%geodAngle))
        firstProfile=firstProfileAsArray(1)

        ! Well, the last profile has to be noInstances later, check this would be OK
        lastProfile=firstProfile+quantity%template%noInstances-1
        if ( lastProfile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Quantity has profiles beyond the end of the l2gp' )

        ! Now check that geodAngle's are a sufficient match
        if ( any(abs(l2gp%geodAngle(firstProfile:lastProfile)-&
          &         quantity%template%phi(1,:)) > tolerance) ) then
          if ( switchDetail(switches, 'l2gp') > -1) then
            call dump ( l2gp%geodAngle(firstProfile:lastProfile), 'L2GP geodetic angle' )
            call dump ( quantity%template%phi(1,:), 'Quantity Geodetic angle' )
          endif
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Quantity has profiles that mismatch l2gp in geodetic angle; interpolate?' )
        end if

        ! Now check that the times match
        if ( any(abs(l2gp%time(firstProfile:lastProfile)- &
          &         quantity%template%time(1,:)) > timeTol) ) then
          if ( switchDetail(switches, 'l2gp') > -1) then
            call dump ( l2gp%time(firstProfile:lastProfile), 'L2GP geodetic angle' )
            call dump ( quantity%template%time(1,:), 'Quantity Geodetic angle' )
          endif
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Quantity has profiles that mismatch l2gp in time' )
        endif
        ! Currently the code cannot interpolate in 3 dimensions, wouldn't
        ! be hard to code up, but no need as yet.
        if ( interpolate .and. quantity%template%noChans /= 1 ) then
          errorCode=cantInterpolate3D
          go to 9
        end if
      else
        ! Given a specific profile, check it's legal
        if ( profile == 0 .or. profile > l2gp%nTimes ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal profile request in l2gp fill' )
      end if

      ! OK, things seem to be OK, so start getting the data
      if ( fromPrecision ) then
        source => l2gp%l2gpPrecision
      else
        source => l2gp%l2gpValue
      end if

      ! OK, now do the filling, it's easier if we don't have to interpolate
      if ( .not. interpolate ) then
        if ( profile == -1 ) then
          ! Not forcing a particular profile to all instances
          quantity%values=reshape(source(:,:,firstProfile:lastProfile),&
            & (/quantity%template%noChans*quantity%template%noSurfs,&
            &   quantity%template%noInstances/))
        else
          ! Spread one profile onto all instances
          quantity%values = spread ( reshape ( source(:,:,profile), &
            &   (/ quantity%template%noChans*quantity%template%noSurfs/) ), &
            & 2, quantity%template%noInstances )
        end if
      else
        if ( quantity%template%verticalCoordinate == l_pressure ) then
          outZeta = -log10 ( quantity%template%surfs(:,1) )
        else
          outZeta = quantity%template%surfs(:,1)
        end if
        ! Setup the interpolation we'll be doing.
        call InterpolateArraySetup ( -log10(real(l2gp%pressures, r8)), &
          & outZeta, 'Linear', coeffs, extrapolate='Clamp' )
        do instance = 1, quantity%template%noInstances
          if ( profile == -1 ) then
            thisProfile = firstProfile + instance - 1
          else
            thisProfile = profile
          end if
          ! Guess I don't really need the loop here, but it
          ! does make the spread stuff much easier, and the setup/teardown
          ! at least makes it more efficient.
          call InterpolateValues ( coeffs, &
            & -log10(real(l2gp%pressures, r8)), &  ! Old X
            & real(source(1,:,thisProfile), r8), & ! OldY
            & outZeta, & ! New X
            & quantity%values(:,instance), & ! New Y
            & method='Linear', extrapolate='Clamp' )
        end do
        call InterpolateArrayTeardown ( coeffs )
      end if
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromL2GP' )

    end subroutine FromL2GP

    ! -------------------------------------- FromProfile --
    subroutine FromProfile_values ( quantity, inheights, invalues, &
      & instancesNode, globalUnit, ptan, logSpace )
      ! Given two arrays,  (heights, invalues) vs. the quantity's own heights,
      ! does the linear interpolation appropriate to perform the fill.
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      real(r8), dimension(:), intent(in) :: inheights
      real(r8), dimension(:), intent(in) :: invalues
      integer, intent(in) :: INSTANCESNODE ! Tree node for instances
      integer, intent(in) :: GLOBALUNIT   ! Possible global unit
      type (VectorValue_T), pointer :: PTAN ! press. values
      logical, intent(in), optional :: LOGSPACE ! Interpolate in logspace

      ! Local variables
      integer :: C                      ! Channel loop counter
      logical, dimension(:), pointer :: DUPLICATED ! Flags
      logical :: Fail                   ! Status from Hunt
      real (r8), dimension(:), pointer :: HEIGHTS ! Heights for the points
      integer :: I,J                    ! Loop counters / indices
      integer, dimension(:), pointer :: ININDS ! Indices
      logical, dimension(:), pointer :: INSTANCES ! Flags
      logical :: LOCALOUTHEIGHTS ! Set if out heights is our own variable
      logical :: MYLOGSPACE             ! Interpolate in log space?
      integer :: NOPOINTS               ! Number of points supplied
      integer :: NOUNIQUE               ! Number of unique heights supplied
      real (r8), dimension(:), pointer :: OUTHEIGHTS ! Heights for output
      real (r8), dimension(:), pointer :: OUTVALUES ! Single profile for output
      integer :: S                      ! Surface loop counter
      integer :: STATUS                 ! Flag
      real (r8), dimension(:), pointer :: VALUES ! Values for the points

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromProfile_values', instancesNode )
      call MLSMessageCalls( 'push', constantName='FromProfile_values' )
      ! Set some stuff up
      myLogSpace = quantity%template%logBasis
      if ( present ( logSpace ) ) myLogSpace = logSpace
      if ( myLogSpace .and. any ( invalues <= 0.0 ) ) then
        call Announce_Error ( 0, no_error_code, &
          & 'Non-positive input data in log profile fill (reset logSpace=false?)' )
        go to 9
      end if
      nullify ( heights, values, duplicated, outHeights, outValues, instances )
      noPoints = size(invalues)
      call Allocate_test ( heights, noPoints, 'values', ModuleName )
      call Allocate_test ( values, noPoints, 'values', ModuleName )
      call Allocate_test ( duplicated, noPoints, 'duplicated', ModuleName )
      call Allocate_test ( outValues, quantity%template%noSurfs, &
        & 'outValues', ModuleName )
      call Allocate_test ( instances, quantity%template%noInstances, &
        & 'instances', ModuleName )
      heights = inHeights
      values = invalues
      if ( myLogSpace ) values = log ( invalues )

      ! Get the appropriate height coordinate for output, for pressure take log.
      if ( associated(ptan) ) then
        localOutHeights = .false.
        outHeights => ptan%values(:,1)
      elseif ( quantity%template%verticalCoordinate == l_pressure ) then
        localOutHeights = .true.
        call Allocate_test ( outHeights, quantity%template%noSurfs, &
          & 'outHeights', ModuleName )
        outHeights = -log10 ( quantity%template%surfs(:,1) )
      else
        localOutHeights = .false.
        outHeights => quantity%template%surfs(:,1)
      end if

      ! Now, if the quantity is coherent, let's assume the user wanted the
      ! 'nearest' values
      if ( quantity%template%coherent .or. associated(ptan) ) then
        nullify ( inInds )
        call allocate_test ( inInds, noPoints, 'inInds', ModuleName )
        ! Hunt fails with non-monotonic outHeights
        if ( .not. isMonotonic(outHeights) ) then
          call monotonize( outHeights )
          if ( WARNWHENPTANNONMONOTONIC ) then                               
            call MLSMessage ( MLSMSG_Warning, ModuleName // '/FromProfile', &
              & 'Ptan non-monotonic' )                                       
            call dump( outHeights, 'outHeights' )
          endif
        endif
        call hunt ( outHeights, heights, inInds, &
         & nearest=.true., allowTopValue=.true., fail=fail )
        if ( fail ) then
          call Announce_Error ( 0, no_error_code, &
          & 'Problem in Hunt' )
          go to 9
        end if
        duplicated = .false.
        do i = 1, noPoints - 1
          do j = i + 1, noPoints
            if ( inInds(i) == inInds(j) ) then
              duplicated ( j ) = .true.
            end if
          end do
        end do
        noUnique = count ( .not. duplicated )
        inInds(1:noUnique) = pack ( inInds, .not. duplicated )
        heights(1:noUnique) = outHeights ( inInds(1:noUnique) )
        values(1:noUnique) = pack ( values, .not. duplicated )
        call deallocate_test ( inInds, 'inInds', ModuleName )
      end if

      ! Now do the interpolation for the first instance
      call InterpolateValues ( heights(1:noUnique), values(1:noUnique), outHeights, &
        & outValues, 'Linear', extrapolate='Constant' )

      if ( myLogSpace ) outValues = exp ( outValues )

      ! Work out what instances we're going to spread this to
      instances = .true.
      if ( instancesNode /= 0 ) then
        call GetIndexFlagsFromList ( instancesNode, instances, &
          & status, noError = .true. )
        if ( status /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to parse instances (wrong units?)' )
      end if

      ! Spread it into the other instances, worry about the mask and instances
      do i = 1, quantity%template%noInstances
        if ( instances(i) ) then
          j = 1
          do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
              if ( associated(quantity%mask) ) then
                if ( iand ( ichar(quantity%mask(j,i)), m_fill ) == 0 ) &
                  & quantity%values(j,i) = outValues(s)
              else
                quantity%values(j,i) = outValues(s)
              end if
              j = j + 1
            end do
          end do
        end if
      end do

      ! Finish off
      if ( localOutHeights ) call Deallocate_test ( outHeights, &
        & 'outHeights', ModuleName )
      call Deallocate_test ( heights, 'heights', ModuleName )
      call Deallocate_test ( values, 'values', ModuleName )
      call Deallocate_test ( outvalues, 'outvalues', ModuleName )
      call Deallocate_test ( duplicated, 'duplicated', ModuleName )
      call Deallocate_test ( instances, 'instances', ModuleName )
      call MLSMessageCalls( 'pop' )
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromProfile_values' )
    end subroutine FromProfile_values

    subroutine FromProfile_node ( quantity, valuesNode, &
      & instancesNode, globalUnit, ptan, logSpace )
      ! This fill is slightly complicated.  Given a values array along
      ! the lines of [ 1000mb : 1.0K, 100mb : 1.0K,  10mb : 2.0K] etc. it
      ! does the linear interpolation appropriate to perform the fill.
      type (VectorValue_T), intent(inout) :: QUANTITY ! Quantity to fill
      integer, intent(in) :: VALUESNODE   ! Tree node for values
      integer, intent(in) :: INSTANCESNODE ! Tree node for instances
      integer, intent(in) :: GLOBALUNIT   ! Possible global unit
      type (VectorValue_T), pointer :: PTAN ! press. values
      logical, intent(in), optional :: LOGSPACE ! Interpolate in logspace

      ! Local variables
      integer :: HEIGHTUNIT             ! Unit for height
      integer :: NOPOINTS               ! Number of points supplied
      integer :: I                      ! Loop counters / indices
      integer :: TESTUNIT               ! Unit for value
      logical :: MYLOGSPACE             ! Interpolate in log space?
      real (r8), dimension(:), pointer :: HEIGHTS ! Heights for the points
      real (r8), dimension(:), pointer :: VALUES ! Values for the points
      real (r8), dimension(2) :: EXPRVALUE ! Value of expression
      integer, dimension(2) :: EXPRUNIT   ! Unit for expression

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.FromProfile_node', instancesNode )
      call MLSMessageCalls( 'push', constantName='FromProfile_node' )

      ! Check the quantity is amenable to this type of fill
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true. ) .and. .not. associated(ptan) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'The quantity is not amenable to a profile fill unless you supply ptan' )

      ! Check the units
      testUnit = quantity%template%unit
      if ( globalUnit /= phyq_Invalid ) testUnit = globalUnit

      ! Set some stuff up
      myLogSpace = quantity%template%logBasis
      if ( present ( logSpace ) ) myLogSpace = logSpace
      noPoints = nsons ( valuesNode ) - 1
      nullify ( heights, values )
      call Allocate_test ( heights, noPoints, 'heights', ModuleName )
      call Allocate_test ( values, noPoints, 'values', ModuleName )

      ! Loop over the values
      do i = 1, noPoints
        ! Get value from tree
        call expr ( subtree ( i+1, valuesNode ), exprUnit, exprValue )
        ! Check height unit OK
        heightUnit = GetUnitForVerticalCoordinate ( quantity%template%verticalCoordinate )
        if ( associated(ptan) ) heightUnit = PHYQ_Zeta
        if ( exprUnit(1) /= heightUnit .and. exprUnit(1) /= PHYQ_Dimensionless &
          & .and. .not. ( exprUnit(1) == PHYQ_Pressure .and. heightUnit == PHYQ_Zeta ) ) &
          & call Announce_error ( valuesNode, no_Error_Code, 'Bad height units for profile fill' )
        ! Store height
        if ( heightUnit == PHYQ_Zeta ) then
          ! Assume zeta coordinates are expressed in mb
          heights(i) = -log10 ( exprValue(1) )
        else
          heights(i) = exprValue(1)
        end if
        ! Check value unit OK
        if ( all ( exprUnit(2) /= (/ testUnit, PHYQ_Dimensionless /) ) ) &
          & call Announce_error ( valuesNode, no_error_code, 'Bad units for profile fill' )
        ! Store value
        values ( i ) = exprValue(2)
      end do

      call FromProfile( quantity, heights, values, &
      & instancesNode, globalUnit, ptan, logSpace )
      call Deallocate_test ( heights, 'heights', ModuleName )
      call Deallocate_test ( values, 'values', ModuleName )
      call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.FromProfile_node' )
    end subroutine FromProfile_node

    ! ---------------------------------------------- ManipulateVectors -----
    subroutine ManipulateVectors ( MANIPULATION, DEST, A, B, C, BOOLEANNAME )
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: GETHASHELEMENT
    use MLSSTRINGS, only: STREQ
      ! Manipulate common items in a, b, copying result to those in dest
      integer, intent(in) :: MANIPULATION
      type (Vector_T), intent(in)            :: A, B
      type (Vector_T), intent(inout)         :: DEST
      real(rv), optional                     :: C  ! constant "c" in manipulation
      character(len=*), intent(in), optional :: BOOLEANNAME

      ! Local variables
      type (VectorValue_T), pointer     :: DQ     ! Destination quantity
      type (VectorValue_T), pointer     :: AQ, BQ ! Source quantities
      integer                           :: SQI    ! Quantity index in source
      integer                           :: N
      character(len=64), dimension(128) :: NAMES
      character(len=64)                 :: QNAME

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ManipulateVectors' )
      n = 0
      ! If we're given a BooleanName, try to interpret it as a container for
      ! quantity names
      if ( present(booleanName) ) then
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & lowerCase(trim(booleanName)), names, countEmpty, &
          & inseparator=runTimeValues%sep )
        n = FindLast( names /= runTimeValues%sep )
      endif
      ! First copy those things in source, loop over them
      dest%globalUnit = a%globalUnit
      do sqi = 1, size ( a%quantities )
        ! Try to find this in dest
        aq => a%quantities(sqi)
        bq => b%quantities(sqi)
        dq => GetVectorQtyByTemplateIndex ( dest, a%template%quantities(sqi) )
        if ( associated(aq) ) call get_string ( aq%template%name, qName, &
          & strip=.true. )
        if ( n > 0 ) then
          if ( .not. any( streq(names, qName, options='-cf' ) )  ) cycle
        elseif ( all( (/ associated(aq), associated(bq), associated(dq) /) ) ) then
          call ByManipulation ( dq, aq, bq, &
            & manipulation, key=0, &
            & ignoreTemplate=.true., spreadflag=.false., dimList=' ', &
            & c=c )
        else
          call MLSMessage ( MLSMSG_Warning, ModuleName // '%ManipulateVectors', &
          & 'dq not associated' )
        end if
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ManipulateVectors' )
    end subroutine ManipulateVectors

    ! -------------------------------------------- WithBinResults -----
    subroutine WithBinResults ( key, quantity, sourceQuantity, ptanQuantity, &
      & channel, method, additional, excludeBelowBottom, centerVertically )
      ! This fills a coherent quantity with the max/min binned value of
      ! a typically incoherent one.  The bins are centered horizontally
      ! on the profiles in quantity, but vertically the bins run between one
      ! surface and the next one up.
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: SOURCEQUANTITY
      type (VectorValue_T), pointer :: PTANQUANTITY
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: METHOD
      logical, intent(in) :: ADDITIONAL
      logical, intent(in) :: EXCLUDEBELOWBOTTOM
      logical, intent(in) :: CENTERVERTICALLY

      ! Local variables
      real(r8), dimension(:,:), pointer :: SOURCEHEIGHTS ! might be ptan.
      real(r8), dimension(:), pointer :: PHIBOUNDARIES
      real(r8) :: V                     ! A value
      integer, dimension(:,:), pointer :: SURFS ! Surface mapping source->quantity
      integer, dimension(:,:), pointer :: INSTS ! Instance mapping source->quantity
      integer :: QS, QI, SS, SI                   ! Loop counters
      integer :: MYCHANNEL              ! Channel or 1
      integer :: NumQtyInstances        ! E.g., num MAFs
      integer :: NumSourceInstances     ! E.g., num MAFs
      logical :: ExtraProfile           ! True when profile stands outside chunk

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithBinResults', key )
      call MLSMessageCalls( 'push', constantName='WithBinResults' )

      ! Check the output quantity
      if ( .not. ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true. ) ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Illegal quantity for bin min/max/total fill' )

      ! Also should have the condition:
      !  quantityType = (/ sourceQuantity%template%quantityType /)
      ! However, the code does not yet really support the ability to do that.
      ! We are alert for the case when HGrid inserted an otherwise
      ! "dropped" profile
      ! This profile will have a phi greater than the chunk-bounding last MAF's phi
      NumSourceInstances = size(sourceQuantity%template%phi, 2)
      NumQtyInstances = size(Quantity%template%phi, 2)
      ExtraProfile = ( quantity%template%noInstances > 1 ) .and. &
        & ( sourceQuantity%template%phi(1,NumSourceInstances) < &
        & quantity%template%phi(1, NumQtyInstances) )

      ! Work out source vertical coordinate
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) then
        nullify ( sourceHeights )
        call Allocate_test ( sourceHeights, sourceQuantity%template%nosurfs, &
          & sourceQuantity%template%noinstances, 'sourceHeights', ModuleName )
        sourceHeights = ptanQuantity%values
        if ( quantity%template%verticalCoordinate /= l_zeta .and. WARNIFVERTCOORDNOTZETA ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinate in quantity to fill is not zeta' )
      else
        sourceHeights => sourceQuantity%template%surfs
        if ( sourceQuantity%template%verticalCoordinate /= quantity%template%verticalCoordinate ) then
          call dump ( quantity%template )
          call dump ( sourceQuantity%template )
          call outputNamedValue( 'source is minor frame?', sourceQuantity%template%minorFrame )
          call outputNamedValue( 'associated(ptanQuantity)?', associated(ptanQuantity) )
          call Announce_Error ( key, no_error_code, &
          & 'Vertical coordinates in binned fill do not match' )
        end if
      end if

      ! Setup index arrays
      nullify ( surfs, insts )
      call Allocate_test ( surfs, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'surfs', ModuleName )
      call Allocate_test ( insts, sourceQuantity%template%nosurfs, &
        & sourceQuantity%template%noinstances, 'insts', ModuleName )

      ! Work out the vertical mapping, a function of instance for
      ! incoherent quantities
      if ( sourceQuantity%template%coherent ) then
        call Hunt ( quantity%template%surfs(:,1), sourceHeights(:,1), surfs(:,1), &
          & allowTopValue=.true., allowBelowValue=excludeBelowBottom, nearest=centerVertically )
        surfs = spread ( surfs(:,1), 2, sourceQuantity%template%noInstances )
      else
        do si = 1, sourceQuantity%template%noInstances
          call Hunt ( quantity%template%surfs(:,1), sourceHeights(:,si), surfs(:,si), &
            & allowTopValue=.true., allowBelowValue=excludeBelowBottom, nearest=centerVertically )
        end do
      end if

!      call dump ( surfs, 'surfs' )

      if ( quantity%template%noInstances > 1 ) then
        ! Setup and work out the phiBoundaries
        nullify ( phiBoundaries )
        call Allocate_test ( phiBoundaries, quantity%template%noInstances-1, &
          & 'phiBoundaries', ModuleName )
        phiBoundaries = 0.5 * ( &
          & quantity%template%phi(1,1:quantity%template%noInstances-1) + &
          & quantity%template%phi(1,2:quantity%template%noInstances) )
        if ( sourceQuantity%template%stacked ) then
          call Hunt ( phiBoundaries, sourceQuantity%template%phi(1,:), &
            & insts(1,:), allowTopValue=.true., allowBelowValue=.true. )
          insts = spread ( surfs(1,:), 1, sourceQuantity%template%noSurfs ) + 1
        else
          do ss = 1, sourceQuantity%template%noSurfs
            call Hunt ( phiBoundaries, sourceQuantity%template%phi(ss,:), &
              & insts(ss,:), allowTopValue=.true., allowBelowValue=.true. )
          end do
          insts = insts + 1
        end if
        if ( method==l_binTotal .and. .false. ) then
          call outputNamedValue( 'phiBoundaries', phiBoundaries )
          call outputNamedValue( 'phi(source)', sourceQuantity%template%phi(1,:) )
        endif
        call Deallocate_test ( phiBoundaries, 'phiBoundaries', ModuleName )
      else
        insts = 1
      end if

!      call dump ( insts, 'insts' )

      ! Modify the surfs index to account for the channel, so now it's more
      ! of a 'values' index
      if ( sourceQuantity%template%frequencyCoordinate /= l_none .and. &
        & channel == 0 ) then
        call Announce_Error ( key, no_error_code, &
          & 'Must supply channel for this bin max/min fill' )
        go to 9
      end if
      myChannel = channel
      if ( channel == 0 ) myChannel = 1

      ! Now loop over the output quantity points and work out the information
      if ( ExtraProfile ) then
        do qs = 1, quantity%template%noSurfs
          if ( count ( surfs == qs ) > 0 ) then
            select case ( method )
            case  ( l_binMax )
              ! Compute the maximum value in the bin
              v = maxval ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, 1 ), &
                & surfs(:,1) == qs ) )
              if ( additional ) then
                quantity%values(qs,:) = max ( quantity%values(qs,:), v )
              else
                quantity%values(qs,:) = v
              end if
            case ( l_binMean )
              ! Compute the average in the bin, be careful about dividing by zero
              v = sum ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, 1 ), &
                & surfs(:,1) == qs ) ) / &
                & max ( count ( surfs(:,1) == qs ), 1 )
              if ( additional ) then
                quantity%values(qs,:) = 0.5 * ( quantity%values(qs,:) + v )
              else
                quantity%values(qs,:) = v
              end if
            case ( l_binMin )
              ! Compute the minimum value in the bin
              v = minval ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, 1 ), &
                & surfs(:,1) == qs ) )
              if ( additional ) then
                quantity%values(qs,:) = min ( quantity%values(qs,:), v )
              else
                quantity%values(qs,:) = v
              end if
            case ( l_binTotal )
              ! Compute the total in the bin
              v = sum ( pack ( sourceQuantity%values ( &
                & myChannel : sourceQuantity%template%instanceLen : &
                &   sourceQuantity%template%noChans, 1 ), &
                & surfs(:,1) == qs ) )
              if ( additional ) then
                quantity%values(qs,:) = quantity%values(qs,:) + v
              else
                quantity%values(qs,:) = v
              end if
            end select
          else
            quantity%values(qs,:) = 0.0
          end if
        end do
      else
        do qi = 1, quantity%template%noInstances
          do qs = 1, quantity%template%noSurfs
            if ( count ( surfs == qs .and. insts == qi ) > 0 ) then
              select case ( method )
              case  ( l_binMax )
                ! Compute the maximum value in the bin
                v = maxval ( pack ( sourceQuantity%values ( &
                  & myChannel : sourceQuantity%template%instanceLen : &
                  &   sourceQuantity%template%noChans, : ), &
                  & surfs == qs .and. insts == qi ) )
                if ( additional ) then
                  quantity%values(qs,qi) = max ( quantity%values(qs,qi), v )
                else
                  quantity%values(qs,qi) = v
                end if
              case ( l_binMean )
                ! Compute the average in the bin, be careful about dividing by zero
                v = sum ( pack ( sourceQuantity%values ( &
                  & myChannel : sourceQuantity%template%instanceLen : &
                  &   sourceQuantity%template%noChans, : ), &
                  & surfs == qs .and. insts == qi ) ) / &
                  & max ( count ( surfs == qs .and. insts == qi ), 1 )
                if ( additional ) then
                  quantity%values(qs,qi) = 0.5 * ( quantity%values(qs,qi) + v )
                else
                  quantity%values(qs,qi) = v
                end if
              case ( l_binMin )
                ! Compute the minimum value in the bin
                v = minval ( pack ( sourceQuantity%values ( &
                  & myChannel : sourceQuantity%template%instanceLen : &
                  &   sourceQuantity%template%noChans, : ), &
                  & surfs == qs .and. insts == qi ) )
                if ( additional ) then
                  quantity%values(qs,qi) = min ( quantity%values(qs,qi), v )
                else
                  quantity%values(qs,qi) = v
                end if
              case ( l_binTotal )
                ! Compute the total in the bin
                v = sum ( pack ( sourceQuantity%values ( &
                  & myChannel : sourceQuantity%template%instanceLen : &
                  &   sourceQuantity%template%noChans, : ), &
                  & surfs == qs .and. insts == qi ) )
                if ( additional ) then
                  quantity%values(qs,qi) = quantity%values(qs,qi) + v
                else
                  quantity%values(qs,qi) = v
                end if
              end select
            else
              quantity%values(qs,qi) = 0.0
            end if
          end do
        end do
      endif
      ! Now tidy up
      call Deallocate_test ( surfs, 'surfs', ModuleName )
      call Deallocate_test ( insts, 'insts', ModuleName )
      if ( associated ( ptanQuantity ) .and. sourceQuantity%template%minorFrame ) &
        & call Deallocate_test ( sourceHeights, 'sourceHeights', ModuleName )
    9 call MLSMessageCalls( 'pop' )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithBinResults' )
    end subroutine WithBinResults

    ! --------------------------------------------- WithBoxcarFunction  ----
    subroutine WithBoxcarFunction ( key, quantity, sourceQuantity, &
    & width, method, ignoreTemplate )
      integer, intent(in) :: KEY        ! Key for tree node
      type (VectorValue_T), intent(inout)      :: QUANTITY
      type (VectorValue_T), intent(in), target :: SOURCEQUANTITY
      integer, intent(in)                      :: WIDTH
      integer, intent(in)                      :: METHOD  ! L_MEAN, L_MAX, L_MIN
      logical, intent(in)              ::         IGNORETEMPLATE
      ! Local variables
      integer :: I, I1, I2              ! Instance indices
      integer :: HALFWIDTH
      real(r8), dimension(:,:), pointer :: OLDVALUES

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.WithBoxcarFunction', key )
      if ( .not. ignoreTemplate .and. &
        & quantity%template%name /= sourceQuantity%template%name ) then
        call Announce_Error ( key, no_error_code, 'Quantity and source quantity do not match' )
        go to 9
      end if
      if ( width <= 1 .or. mod ( width, 2 ) == 0 ) then
        call Announce_Error ( key, no_error_code, 'width must be greater than 1 and odd' )
        go to 9
      end if

      if ( associated ( quantity%values, sourceQuantity%values ) ) then
        nullify ( oldValues )
        call Allocate_test ( oldValues, quantity%template%instanceLen, quantity%template%noInstances, &
          & 'oldValues', ModuleName )
        oldValues = sourceQuantity%values
      else
        oldValues => sourceQuantity%values
      end if

      halfWidth = width/2
      select case ( method )
      case ( l_mean )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = sum ( oldValues(:,i1:i2), dim=2 ) / (i2-i1+1)
        end do
      case ( l_min )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = minval ( oldValues(:,i1:i2), dim=2 )
        end do
      case ( l_max )
        do i = 1, quantity%template%noInstances
          i1 = max ( i - halfWidth, 1 )
          i2 = min ( i + halfWidth, quantity%template%noInstances )
          quantity%values(:,i) = maxval ( oldValues(:,i1:i2), dim=2 )
        end do
      end select

      if ( associated ( quantity%values, sourceQuantity%values ) ) then
        call Deallocate_test ( oldValues, 'oldValues', ModuleName )
      end if
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.WithBoxcarFunction' )

    end subroutine WithBoxcarFunction

    ! ----------------------------------------- offsetradiancequantity -----
    subroutine OffsetRadianceQuantity ( quantity, radianceQuantity, amount )
      type (VectorValue_T), intent(inout) :: QUANTITY
      type (VectorValue_T), intent(in) :: RADIANCEQUANTITY
      real (rv), intent(in) :: AMOUNT

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.OffsetRadianceQuantity' )
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity for offsetRadiance fill is not radiance' )
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Radiance quantity for offsetRadiance fill is not radiance' )
      if ( quantity%template%signal /= radianceQuantity%template%signal .or. &
        & quantity%template%sideband /= radianceQuantity%template%sideband ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity and rad. qty. in offsetRadiance fill different signal/sideband' )
      if ( associated ( radianceQuantity%mask ) ) then
        where ( iand ( ichar(radianceQuantity%mask), M_LinAlg ) /= 0 )
          quantity%values = quantity%values + amount
        end where
      end if
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.OffsetRadianceQuantity' )
    end subroutine OffsetRadianceQuantity

    ! ---------------------------------------------- ResetUnusedRadiances --
    subroutine ResetUnusedRadiances ( quantity, amount )
      type (VectorValue_T), intent(inout) :: QUANTITY
      real (rv), intent(in) :: AMOUNT
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.ResetUnusedRadiances' )
      if ( .not. ValidateVectorQuantity ( quantity, &
        & quantityType=(/l_radiance/) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Quantity for resetUnusedRadiances fill is not radiance' )
      where ( quantity%values > amount*0.9 )
        quantity%values = quantity%values - amount
      end where
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.ResetUnusedRadiances' )
    end subroutine ResetUnusedRadiances

    ! --------------------------------------------- RotateMagneticField ----
    subroutine RotateMagneticField ( key, qty, fieldECR, ecrToFOV )
      integer, intent(in) :: KEY        ! Where are we in the l2cf?
      type (VectorValue_T), intent(inout) :: QTY ! The quantity to fill
      type (VectorValue_T), intent(in) :: FIELDECR ! The input field
      type (VectorValue_T), intent(in) :: ECRTOFOV ! The rotation matrix
      ! Local variables
      integer :: INSTANCE               ! Loop counter
      integer :: SURFACE                ! Loop counter
      integer :: MAF(1)                 ! Which MAF is the best match to this instance
      real(r8) :: THISFIELD(3)          ! A magnetic field vector
      real(r8) :: ROTATION(3,3)         ! One rotation matrix
      real(r8) :: STRENGTH              ! The field strength

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.RotateMagneticField', key )
      ! Do some sanity checks
      if ( .not. any ( qty%template%quantityType == &
        & (/ l_fieldStrength, l_fieldAzimuth, l_fieldElevation /) ) ) then
        call Announce_Error ( key, no_error_code, 'Inappropriate quantity for this fill' )
        go to 9
      end if
      if ( .not. DoHGridsMatch ( qty, fieldECR ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Field and result quantity must have matching hGrids' )
        go to 9
      end if
      if ( .not. DoVGridsMatch ( qty, fieldECR ) ) then
        call Announce_Error ( key, no_error_code, &
          & 'Field and result quantity must have matching hGrids' )
        go to 9
      end if

      do instance = 1, qty%template%noInstances
        ! Identify the relevant MAF and pull out the first MIF's rotation matrix
        ! (first MIF is probably close enough to all the others to be useful).
        maf = minloc ( abs ( qty%template%phi(1,instance) - &
          & ecrToFOV%template%phi(1,:) ) )
        rotation = reshape ( ecrToFOV%values( 1:9, maf), (/ 3, 3 /) )
        ! Now loop over the pressure levels in the input and output field information
        do surface = 1, qty%template%noSurfs
          ! Get the field in ECR coordinates
          thisField = fieldECR%values ( (surface-1) * 3 + 1 : surface*3, &
            & instance )
          ! Now rotate it into IFOVPP coordinates
          thisField = matmul ( rotation, thisField )
          ! Now work out the strength / angles as appropriate.
          strength = sqrt ( sum ( thisField ** 2 ) )
          select case ( qty%template%quantityType )
          case ( l_fieldStrength )
            qty%values(surface,instance) = strength
          case ( l_fieldElevation )
            if ( strength /= 0.0_r8 ) then
              ! Elevation is constrained to 0--90 degrees instead of 0--180 degrees because
              ! radiative transfer Physics is symmetric.  We save half of l2pc bins.
              qty%values(surface,instance) = acos ( abs ( thisField(3) / strength ) ) * rad2deg
            else
              qty%values(surface,instance) = 0.0
            end if
          case ( l_fieldAzimuth )
            if ( thisField(1) /= 0.0_r8 ) then
              qty%values(surface,instance) = atan ( thisField(2) / thisField(1) ) * rad2deg
            else
              qty%values(surface,instance) = merge ( 90.0_r8, -90.0_r8, &
                & thisField(1) > 0.0_r8 )
            end if
          end select
        end do
      end do
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.RotateMagneticField' )

    end subroutine RotateMagneticField

    ! ---------------------------------------------- ScaleOverlaps -------
    subroutine ScaleOverlaps ( quantity, multiplierNode )
      type (VectorValue_T), intent(inout) :: QUANTITY
      integer, intent(in) :: MULTIPLIERNODE ! Tree node for factors
      ! Local variables
      real (r8) :: exprValue(2)         ! Tree expression
      integer :: exprUnit(2)            ! Tree units
      integer :: i,j                    ! Loop counters / indices

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.ScaleOverlaps', multiplierNode )
      scaleLowerLoop: do i = 1, quantity%template%noInstancesLowerOverlap
        if ( i+1 > nsons ( multiplierNode ) ) exit scaleLowerLoop
        call expr_check ( subtree( i+1, multiplierNode ), exprUnit, exprValue, &
          & (/PHYQ_Dimensionless/), unitsError )
        if ( unitsError ) then
          call Announce_error ( multiplierNode, wrongUnits, &
            & extraMessage="for scaleOverlaps fill", &
            & extraInfo=(/exprUnit(1), PHYQ_Dimensionless/) )
          go to 9
        end if
        if ( associated ( quantity%mask ) ) then
          where ( iand ( ichar(quantity%mask(:,i)), m_Fill ) == 0 )
            quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
          end where
        else
          quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
        end if
      end do scaleLowerLoop

      scaleUpperLoop: do i = quantity%template%noInstances, &
        & quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap + 1, - 1
        j = quantity%template%noInstances - i + 1
        if ( j+1 > nsons ( multiplierNode ) ) exit scaleUpperLoop
        call expr_check ( subtree( j+1, multiplierNode ), exprUnit, exprValue, &
          & (/PHYQ_Dimensionless/), unitsError )
        if ( unitsError ) then
          call Announce_error ( multiplierNode, wrongUnits, &
            & extraMessage="for scaleOverlaps fill", &
            & extraInfo=(/exprUnit(1), PHYQ_Dimensionless/) )
          go to 9
        end if
        if ( associated ( quantity%mask ) ) then
          where ( iand ( ichar(quantity%mask(:,i)), m_Fill ) == 0 )
            quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
          end where
        else
          quantity%values ( :, i ) = quantity%values ( :, i ) * exprValue(1)
        end if
      end do scaleUpperLoop
    9 if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.ScaleOverlaps' )
    end subroutine ScaleOverlaps

    ! ----------------------------------------- SpreadChannelFill --------
    subroutine SpreadChannelFill ( quantity, channel, key, &
      & sourceQuantity )
      type(VectorValue_T), intent(inout) :: QUANTITY
      integer, intent(in) :: CHANNEL
      integer, intent(in) :: KEY
      type(VectorValue_T), intent(in), optional :: SOURCEQUANTITY
      ! Local variables
      integer :: I                      ! Instance loop counter
      integer :: C                      ! Channel loop counter
      integer :: S                      ! Surface loop counter
      integer :: J                      ! Destination index
      integer :: MYCHANNEL              ! Possibly offset channel
      type (Signal_T) ::  signal        ! Signal for this quantity

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.SpreadChannelFill', key )
      if ( present(sourceQuantity) ) then
        do i = 1, quantity%template%noInstances
          do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
              j = (s-1)*quantity%template%noChans + c
              if ( associated ( quantity%mask ) ) then
                if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) == 1 ) cycle
              end if
              quantity%values ( j, i ) = &
                & sourceQuantity%values ( s, i )
            end do
          end do
        end do
      else
        ! Deal with any channel numbering issues.
        signal = GetSignal ( quantity%template%signal )
        myChannel = channel - lbound ( signal%frequencies, 1 ) + 1
        if ( myChannel < 1 .or. myChannel > quantity%template%noChans ) &
          & call Announce_Error ( key, no_error_code, &
          & 'Invalid channel for spread channel' )
        do i = 1, quantity%template%noInstances
          do s = 1, quantity%template%noSurfs
            do c = 1, quantity%template%noChans
              if ( c == myChannel ) cycle
              j = (s-1)*quantity%template%noChans + c
              if ( associated ( quantity%mask ) ) then
                if ( iand ( ichar(quantity%mask(j,i)), m_Fill ) == 1 ) cycle
              end if
              quantity%values ( j, i ) = &
                & quantity%values ( (s-1)*quantity%template%noChans + myChannel, i )
            end do
          end do
        end do
      end if
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.SpreadChannelFill' )
    end subroutine SpreadChannelFill

    ! ---------------------------------------------- TransferVectors -----
    subroutine TransferVectors ( SOURCE, DEST, SKIPMASK, INTERPOLATE, &
      & BOOLEANNAME )
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: GETHASHELEMENT, SWITCHDETAIL
    use MLSSTRINGS, only: STREQ
    use OUTPUT_M, only: OUTPUT, OUTPUTNAMEDVALUE
    use TOGGLES, only: SWITCHES
      ! Copy common items in source to those in dest
      type (Vector_T), intent(in) :: SOURCE
      type (Vector_T), intent(inout) :: DEST
      logical, intent(in) :: SKIPMASK
      logical, intent(in) :: INTERPOLATE
      character(len=*), intent(in), optional :: BOOLEANNAME

      ! Local variables
      type (VectorValue_T), pointer :: DQ ! Destination quantity
      type (VectorValue_T), pointer :: SQ ! Source quantity
      integer                           :: N
      character(len=64), dimension(128) :: NAMES
      character(len=64)                 :: QNAME
      integer :: SQI                      ! Quantity index in source
      logical :: verbose, verboser

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.TransferVectors' )
      verbose = ( switchDetail(switches, 'bool') > -1 )
      verboser = ( switchDetail(switches, 'bool') > 0 )
      n = 0
      ! If we're given a BooleanName, try to interpret it as a container for
      ! quantity names
      if ( present(booleanName) ) then
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & lowerCase(trim(booleanName)), names, countEmpty, &
          & inseparator=runTimeValues%sep )
        n = FindLast( names /= runTimeValues%sep )
        if ( verboser ) then
          call outputNamedValue( 'n', n )
          call dump( names(1:n), 'names', width=1 )
        endif
      endif

      ! First copy those things in source, loop over them
      dest%globalUnit = source%globalUnit
      do sqi = 1, size ( source%quantities )
        ! Try to find this in dest
        sq => source%quantities(sqi)
        dq => GetVectorQtyByTemplateIndex ( dest, source%template%quantities(sqi) )
        if ( .not. associated(dq) ) cycle
        call get_string ( dq%template%name, qName, &
          & strip=.true. )
        if ( n > 0 ) then
          if ( .not. any( streq(names, qName, options='-cf' ) )  ) cycle
        endif
        call output( 'Transferring quantities named ' // trim(qName), advance='yes' )
        dq%values = sq%values
        if ( .not. skipMask ) then
          if ( associated(sq%mask) ) then
            if ( .not. associated(dq%mask)) call CreateMask ( dq )
            dq%mask = sq%mask
          else
            call destroyVectorQuantityMask ( dq )
          end if
        elseif ( interpolate ) then
          dq => GetVectorQuantityByType ( dest, &
            & quantityType=sq%template%quantityType, &
            & molecule=sq%template%molecule )
          if ( associated ( dq ) ) then
            call FromInterpolatedQty( dq, sq, key=0, &
              & ignoreTemplate=.false. )
          end if
        end if
      end do
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.TransferVectors' )
    end subroutine TransferVectors

    ! ---------------------------------------------- TransferVectorsByMethod -----
    ! Fill items in dest according to method
    ! Using either those in source or those in measVector

    ! Note:
    ! This may fill multiple quantities in destination using the same
    ! source quantities if more than one match according to
    ! method-dependent criteria (should we check? warn?)
    
    ! Which methods to permit is a question probably to be revisited
    ! as needed
    subroutine TransferVectorsByMethod ( KEY, DEST, &
      & SOURCE, METHOD, DONTMASK, INTERPOLATE, &
      & IGNORENEGATIVE, IGNOREZERO, MEASVECTOR, MODELVECTOR, &
      & NOISEVECTOR, PTAN, BOOLEANNAME )
    use MLSL2OPTIONS, only: RUNTIMEVALUES
    use MLSSTRINGLISTS, only: GETHASHELEMENT
    use MLSSTRINGS, only: STREQ
      integer, intent(in)            :: KEY
      type (Vector_T), pointer       :: DEST
      type (Vector_T), pointer       :: SOURCE
      integer, intent(in)            :: METHOD
      logical, intent(in)            :: DONTMASK
      logical, intent(in)            :: INTERPOLATE
      logical, intent(in)            :: IGNORENEGATIVE
      logical, intent(in)            :: IGNOREZERO
      type (Vector_T), pointer       :: MODELVECTOR
      type (Vector_T), pointer       :: MEASVECTOR
      type (Vector_T), pointer       :: NOISEVECTOR
      type (VectorValue_T), pointer  :: PTAN ! tangent pressure
      character(len=*), intent(in), optional :: BOOLEANNAME

      ! Local variables
      logical, parameter            :: ASPERCENTAGE = .false.
      integer                       :: BININDEXER
      integer, parameter            :: CHANNEL=0
      type (VectorValue_T), pointer :: DQ ! Destination quantity
      integer                       :: DQI
      logical, parameter            :: IGNORETEMPLATE = .false.
      real, dimension(2), parameter :: MULTIPLIER = (/ 1.0, 1.0 /)
      integer                       :: NUMQUANTITIES
      type (VectorValue_T), pointer :: SQ ! Source or meas quantity
      integer                       :: SQI                      ! Quantity index in source or meas
      type (VectorValue_T), pointer :: MQ ! Model quantity
      integer                           :: N
      character(len=64), dimension(128) :: NAMES
      character(len=64)                 :: QNAME
      type (VectorValue_T), pointer :: NQ ! Noise quantity
      integer                       :: NUMMATCHES
      real(r8)                      :: SCALEINSTANCES, SCALERATIO, SCALESURFS
      integer, dimension(4), parameter :: BINMETHODS = &
        & (/ l_binMax, l_binMean, l_binMin, l_binTotal /)
      integer, dimension(4), parameter :: BINNEDTYPES = &
        & (/ l_chisqBinned, l_chisqBinned, l_cloudMinMax, l_noRadsBinned /)
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.TransferVectorsByMethod' )
      n = 0
      ! If we're given a BooleanName, try to interpret it as a container for
      ! quantity names
      if ( present(booleanName) ) then
        call GetHashElement( runTimeValues%lkeys, runTimeValues%lvalues, &
          & lowerCase(trim(booleanName)), names, countEmpty, &
          & inseparator=runTimeValues%sep )
        n = FindLast( names /= runTimeValues%sep )
      endif

      scaleInstances=-1.0d0
      scaleRatio=1.0d0
      scaleSurfs=-1.0d0
      ! First copy those things in source, loop over them
      if ( any( method == &
        & (/ l_chisqMMAF, l_chisqMMIF, l_chisqChan, l_noRadsPerMif /) ) ) then
        NumQuantities = size ( measVector%quantities )
      else
        NumQuantities = size ( source%quantities )
      endif
      do sqi = 1, NumQuantities
        if ( any( method == &
          & (/ l_chisqMMAF, l_chisqMMIF, l_chisqChan /) ) ) then
          sq => measVector%quantities(sqi)
          mq => modelVector%quantities(sqi)
          nq => noiseVector%quantities(sqi)
        elseif ( any( method == &
          & (/ l_noRadsPerMif /) ) ) then
          sq => measVector%quantities(sqi)
        else
          sq => source%quantities(sqi)
        endif
        nummatches = 0
        do dqi = 1, size ( dest%quantities ) 
          dq => dest%quantities(dqi)
          if ( .not. associated(dq) ) cycle
          call get_string ( dq%template%name, qName, &
            & strip=.true. )
          if ( n > 0 ) then
            if ( .not. any( streq(names, qName, options='-cf' ) )  ) cycle
          endif
          select case ( method )
          case ( l_addNoise ) ! ----- Add random noise to source Quantity -------
            nq => GetVectorQtyByTemplateIndex ( noiseVector, &
              & source%template%quantities(sqi) )
            call addGaussianNoise ( key, dq, sq, &
              & nq, multiplier, spread=.false., &
              & ignoreTemplate=ignoreTemplate )
          case ( l_binMax, l_binMean, l_binMin, l_binTotal )
            if ( dq%template%signal /= sq%template%signal ) cycle
            binIndexer = FindFirst( binMethods, method )
            if ( dq%template%quantityType /= binnedTypes(binIndexer) ) cycle
            if ( .not. any( sq%template%quantityType == &
              & (/ l_chiSQMMIF, l_noRadsPerMif/) ) ) cycle
            call WithBinResults ( key, dq, sq, ptan, &
              & channel, method, additional=.false., &
              & excludeBelowBottom=.false., centerVertically=.false. )
          case ( l_lsLocal, l_lsGlobal, l_lsWeighted )
            call UsingLeastSquares ( key, dq, sq, ptan, &
              & channel, method, &
              & scaleInstances, scaleRatio, scaleSurfs )
          case ( l_boxcar )
            if ( dq%template%quantityType /= l_chisqBinned ) cycle
            call WithBoxcarFunction ( key, dq, sq, width=3, &
              & method=l_min, ignoreTemplate=ignoreTemplate )
          case ( l_chiSQChan )
            if ( dq%template%signal /= sq%template%signal ) cycle
            if ( dq%template%quantityType /= l_chiSQChan ) cycle
            call ChiSqChan ( key, dq, &
              & sq, mq, nq, &
              & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
          case ( l_chiSQMMaf )
            if ( dq%template%signal /= sq%template%signal ) cycle
            if ( dq%template%quantityType /= l_chiSQMMaf ) cycle
            call ChiSqMMaf ( key, dq, &
              & sq, mq, nq, &
              & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
          case ( l_chiSQMMIF )
            if ( dq%template%signal /= sq%template%signal ) cycle
            if ( dq%template%quantityType /= l_chiSQMMIF ) cycle
            call ChiSqMMiF ( key, dq, &
              & sq, mq, nq, &
              & dontMask, ignoreZero, ignoreNegative, ignoreTemplate, multiplier )
          case ( l_noRadsPerMif )
            if ( dq%template%signal /= sq%template%signal ) cycle
            if ( dq%template%quantityType /= l_noRadsPerMif ) cycle
            call NoRadsPerMif ( key, dq, sq, asPercentage )
          case default
           call MLSMessage ( MLSMSG_Error, ModuleName, &
             & "Unable to transfer by this method")
          end select
          nummatches = nummatches + 1
        end do
        if ( nummatches > 1 ) then
          L2CFNode = key
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'More than one matching destination quantity in ' // &
          & 'TransferVectorsByMethod' )
        endif
      end do
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.TransferVectorsByMethod' )
    end subroutine TransferVectorsByMethod

    ! ---------------------------------------------- UncompressRadiance -----
    subroutine UncompressRadiance ( key, quantity, totalPowerQuantity, &
      & systemTemperatureQuantity, termsNode )
      integer, intent(in) :: KEY        ! Tree node for messages
      type (VectorValue_T), intent(inout) :: QUANTITY ! The quantity to mess with
      type (VectorValue_T), intent(inout) :: TOTALPOWERQUANTITY ! The total power
      type (VectorValue_T), intent(inout) :: SYSTEMTEMPERATUREQUANTITY ! The system temperature
      integer, intent(in) :: TERMSNODE

      ! Local parameters
      integer, parameter :: NOTERMS=3

      ! Local variables
      integer :: I
      integer :: J 
      integer :: K
      integer :: Nchannels
      integer :: Npointings
      real(rv), dimension(noTerms) :: myTerms
      integer, DIMENSION(2) :: UNITASARRAY ! Unit for value given
      real (r8), DIMENSION(2) :: VALUEASARRAY ! Value give
      real(rv) :: b
      real(rv) :: bTsys
      real(rv) :: sumCal
      real(rv) :: bprodCal

      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.UncompressRadiance', key )
      if ( quantity%template%quantityType /= l_radiance ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Inappropriate quantity for uncompressRadiance fill' )
      if ( totalPowerQuantity%template%quantityType /= l_baseline  ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Inappropriate quantity for total power quantity in uncompressRadiance' // &
        & ' (should be a baseline)' )
      if ( .not. totalPowerQuantity%template%minorFrame ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Inappropriate quantity for total power quantity in uncompressRadiance' // &
        & ' (should be MIF baseline)' )

     ! if (totalPowerQuantity%template%radiometer /= quantity%template%radiometer ) &
     !   & call Announce_Error ( key, no_error_code, &
     !   & 'Inappropriate quantity for total power quantity in uncompressRadiance' // &
     !   & ' (should be the right radiometer)' )

      if ( nsons ( termsNode ) /= noTerms + 1 ) &
        & call Announce_Error ( key, no_error_code, &
        & 'Wrong number of terms for uncompressRadiance fill' )
      
      do i = 1, noTerms
        call expr_check ( subtree(i+1,termsNode), unitAsArray, valueAsArray, &
          & (/PHYQ_Dimensionless/), unitsError )
        if ( unitsError ) call Announce_error ( termsNode, wrongUnits, &
          & extraInfo=(/unitAsArray(1), PHYQ_Temperature/) )
        myTerms(i) = valueAsArray(1)
      end do

      ! Radiances are in quantity%values [ chans * pointings, mafs ], do your messing here
      ! Total power are in totalPowerQuantity%values [ pointings, mafs ]
      ! System temperature is in systemTemperatureQuantity%values(chans,1)
      ! The terms supplied are in myTerms(:)

      ! YOUR CODE HERE  
      !
      !terms needed:
      !   myTerms(1) == b == gain compression at the cal target, e.g. 0.01 
      !   myTerms(2) == Tcal == cal target Planckified radiance, e.g. 250 
      !   myTerms(3) == Tspace == space Planckified radiance,    e.g. 3
      !
      ! TA ~= TAhat + b * TAhat *(TAhatbar - Tcal - Tspace - Tsys) +  b * Tcal * Tspace  +  b * Tsys * TAhatbar
      !
      Nchannels = size(quantity%values, 1)/ size(totalPowerQuantity%values, 1)
      !Nchannels = quantity% 
      Npointings = size(totalPowerQuantity%values, 1)
      do i = 1, size(quantity%values, 2)  !MAFS
        
        do j = 1, Npointings !Pointings
          if ( totalPowerQuantity%values(j,i) == 0 ) then
             call Announce_Error ( key, noTotalPower )
          else
            do k = 1, Nchannels  !channels
              b=myTerms(1)/(myTerms(2) +  systemTemperatureQuantity%values(k,1))
              sumCal = myTerms(2) + myTerms(3) + systemTemperatureQuantity%values(k,1)
              bprodCal = b * myTerms(2) * myTerms(3) 
              bTsys= b * systemTemperatureQuantity%values(k,1)

              quantity%values(k + (j - 1)*Nchannels, i) = quantity%values(k + (j - 1)*Nchannels, i) &
                       &  + b *  quantity%values(k + (j - 1)*Nchannels, i) * (totalPowerQuantity%values(j,i)-sumcal) &
                       &  + bprodCal  + bTsys * totalPowerQuantity%values(j,i)  
            end do
          end if
        end do
      end do
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.UncompressRadiance' )

    end subroutine UncompressRadiance

    ! ------------------------------------------- QtyFromFile ----------
    subroutine QtyFromFile ( key, quantity, MLSFile, &
      & filetype, options, sdname, spread, interpolate )
      use MLSHDF5, only: MATCHHDF5ATTRIBUTES
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Radiance quantity to modify
      type (MLSFile_T), pointer   :: MLSFile
      character(len=*), intent(in) :: FILETYPE
      character(len=*), intent(in) :: OPTIONS
      character(len=*), optional, intent(in) :: sdname
      logical, intent(in)                    :: spread
      logical, intent(in)                    :: interpolate
      ! Local variables
      integer, parameter                      :: MAXLISTLENGTH=256
      character (LEN=10*MAXLISTLENGTH)        :: attrnames
      character (LEN=10*MAXLISTLENGTH)        :: attrvalues
      logical homogeneous
      character(len=80) :: name
      character (len=80) :: Str
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.QtyFromFile', key )
      homogeneous = index(lowercase(options), 'h') > 0
      ! How do we access the dataset to read? By name or by attribute?
      if ( index(lowercase(options), 'a') < 1 ) then
        ! By name
        name = ' '
        if ( present(sdName) ) name = sdname
        if ( len_trim(name) < 1 ) &
          & call get_string( lit_indices(quantity%template%name), Name )
        if ( addSlash ) Name = '/' // Name
      else
        ! By attribute
        attrnames  = 'TemplateName,tempQtyType'
        call Get_String ( quantity%template%name, str, strip=.true. )
        attrvalues = str
        call Get_String ( lit_indices(quantity%template%quantityType), &
          & str, strip=.true. )
        attrvalues = trim(attrvalues) // ',' // str
        call MatchHDF5Attributes ( MLSFile, attrNames, attrValues, name )
        ! call Announce_Error ( key, no_Error_Code, &
        ! &   'Unable to read Quantity from File by attribute yet' )
      endif
      call NamedQtyFromFile ( key, quantity, MLSFile, &
        & filetype, name, spread, interpolate, homogeneous )
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.QtyFromFile' )
    end subroutine QtyFromFile
    
    ! ------------------------------------------- VectorFromFile ----------
    subroutine VectorFromFile ( key, Vector, MLSFile, &
      & filetype, options, spread, interpolate )
      use MLSHDF5, only: GETALLHDF5DSNAMES, MATCHHDF5ATTRIBUTES
      integer, intent(in) :: KEY        ! Tree node
      type (Vector_T), intent(inout) :: Vector
      type (MLSFile_T), pointer   :: MLSFile
      character(len=*), intent(in) :: FILETYPE
      character(len=*), intent(in) :: OPTIONS
      logical, intent(in)                    :: spread
      logical, intent(in)                    :: interpolate
      ! Local variables
      integer :: DSI                      ! Dataset index in file
      character(len=MAXSDNAMESBUFSIZE) :: DSNames
      logical :: homogeneous
      character(len=64) :: name
      type (VectorValue_T), pointer :: quantity 
      integer, parameter                      :: MAXLISTLENGTH=256
      character (LEN=10*MAXLISTLENGTH)        :: attrnames
      character (LEN=10*MAXLISTLENGTH)        :: attrvalues
      integer :: SQI                      ! Quantity index in source
      character (len=80) :: Str
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_begin ( 'FillUtils_1.VectorFromFile', key )
      homogeneous = index(lowercase(options), 'h') > 0
      call output( 'Now in VectorFromFile', advance='yes' )
      call GetAllHDF5DSNames( MLSFile, DSNames )
      call dump( DSNames )
      do dsi=1, NumStringElements( trim(DSNames), countEmpty )
        do sqi = 1, size ( vector%quantities )
          quantity => vector%quantities(sqi)
          ! How do we access the dataset to read? By name or by attribute?
          if ( index(lowercase(options), 'a') < 1 ) then
            ! By name
            if ( index(lowercase(options), 'num') < 1 ) then
              if ( quantity%template%name < 1 ) &
                & call Announce_Error ( key, no_Error_Code, &
                &   'template name is 0?' )
              call get_string( quantity%template%name, Name )
            else
              call writeIntsToChars ( sqi, Name )
              Name = 'Quantity ' // trim(Name)
            endif
            if ( lowercase(trim(name)) /= &
              & lowercase(StringElement( DSNames, dsi, countEmpty )) ) &
              & cycle
            if ( addSlash ) Name = '/' // Name
            if ( DEEBUG ) then
              call outputNamedValue( 'shape(values)' // trim(name), &
                & shape(quantity%values) )
              call dump( MLSFile )
            endif
          else
            ! By attribute
            attrnames  = 'TemplateName,tempQtyType'
            call Get_String ( quantity%template%name, str, strip=.true. )
            attrvalues = str
            call Get_String ( lit_indices(quantity%template%quantityType), &
              & str, strip=.true. )
            attrvalues = trim(attrvalues) // ',' // str
            call MatchHDF5Attributes ( MLSFile, attrNames, attrValues, name )
            ! call Announce_Error ( key, no_Error_Code, &
            ! &   'Unable to read Vector from File by attribute yet' )
          endif
          if ( len_trim(name) > 0 ) &
            & call NamedQtyFromFile ( key, quantity, MLSFile, &
            & filetype, name, spread, interpolate, homogeneous )
        enddo
      enddo
      if ( toggle(gen) .and. levels(gen) > 1 ) &
        & call trace_end ( 'FillUtils_1.VectorFromFile' )
    end subroutine VectorFromFile

  ! -- Private procedures ----------

    ! ------------------------------------------- FillableChiSq ---
    function FillableChiSq ( qty, measQty, modelQty, noiseQty ) result ( aok )
      ! Purpose (A)
      !   Check whether we may proceed with special fill of chi squared
      !   case where all VectorValue_T args present
      ! Purpose (B)
      !   Check whether we may proceed with special fill of addNoise
      !   case where missing noiseQty arg
      type (VectorValue_T), intent(in) ::              QTY
      type (VectorValue_T), intent(in) ::              modelQty
      type (VectorValue_T), intent(in) ::              measQty
      type (VectorValue_T), optional, intent(in) ::    noiseQty
      logical ::                                       AOK

      ! What we will check is that (for the args we have been given):
      ! (0) all quantities have associated values
      ! (1) all quantities have same molecule
      ! (2) all quantities have same signal
      ! (3) all quantities have same HGrid
      ! (4A) all but qty (chiSq) have same VGrid
      ! (4B) all have same VGrid
      ! Ah, but radiances have no VGrid or HGrid,
      ! and molecule should be part of signal, so
      ! if radiances need only check on signal
      ! if vmr, check on others

     ! Local variables
      logical ::       minorFrame   ! TRUE if radiances, FALSE if vmr

      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.FillableChiSq' )
      aok = .true.

      ! (0)
      if ( present(noiseQty) ) then
        aok = associated(noiseQty%values)
        if ( DEEBUG .and. .not. associated(noiseQty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Noise values unassociated in FillableChiSq')
      end if
      aok = aok .and. &
        & associated(qty%values) .and. &
        & associated(measQty%values) .and. &
        & associated(modelQty%values)

      if ( DEEBUG ) then
        if ( .not. associated(qty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Quantity values unassociated in FillableChiSq')
        if ( .not. associated(measQty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Measurements values unassociated in FillableChiSq')
        if ( .not. associated(modelQty%values) ) &
          & call announce_error( no_Error_Code, No_Error_code, &
          & 'Model values unassociated in FillableChiSq')
      end if

      if ( .not. aok ) go to 9

      minorFrame = qty%template%minorFrame .or. qty%template%majorFrame
      ! (1)
      if ( .not. minorFrame ) then
        aok = aok .and. &
          & (qty%template%molecule == measQty%template%molecule) &
          & .and. &
          & (qty%template%molecule == modelQty%template%molecule)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%molecule == noiseQty%template%molecule)
      end if

      ! (2)
      if ( minorFrame ) then
        aok = aok .and. &
          & (qty%template%signal == measQty%template%signal) &
          & .and. &
          & (qty%template%signal == modelQty%template%signal)
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & (qty%template%signal == noiseQty%template%signal)
      end if

      ! (3)
      if ( .not. minorFrame ) then
        aok = aok .and. &
          & DoHgridsMatch( qty, measQty ) &
          & .and. &
          & DoHgridsMatch( qty, modelQty )
        if ( present(noiseQty) ) aok = aok &
          & .and. &
          & DoHgridsMatch( qty, noiseQty )

        ! (4)
        aok = aok .and. &
          & DoVgridsMatch( measqty, modelQty )
        if ( present(noiseQty) ) then
          aok = aok &
          & .and. &
          & DoVgridsMatch( measqty, noiseQty )
        else
          aok = aok &
          & .and. &
          & DoVgridsMatch( measqty, Qty )
        end if
      end if

    9 if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.FillableChiSq' )
    end function FillableChiSq

    ! ------------------------------------------- NamedQtyFromFile ----------
    subroutine NamedQtyFromFile ( key, quantity, MLSFile, &
      & filetype, name, spread, interpolate, homogeneous )
      use MLSHDF5, only: GETHDF5ATTRIBUTE, GETHDF5DSDIMS, LOADFROMHDF5DS
      integer, intent(in) :: KEY        ! Tree node
      type (VectorValue_T), intent(inout) :: QUANTITY ! Radiance quantity to modify
      type (MLSFile_T), pointer   :: MLSFile
      character(len=*), intent(in) :: FILETYPE
      character(len=*), intent(in) :: name
      logical, intent(in)                    :: spread
      logical, intent(in)                    :: interpolate
      logical, intent(in)                    :: homogeneous
      ! Local variables
      integer, dimension(3) :: dimInts
      integer(kind=hSize_t), dimension(3) :: dims
      integer, parameter :: globalUnit = 0
      integer, parameter :: instancesNode = 0
      integer :: instance
      type( L2AUXData_T ) :: l2aux ! Result
      type( l2GPData_T ) :: l2gp ! Result
      integer :: noSurfs
      type (VectorValue_T), pointer :: PTAN => null()
      integer :: status
      real(r8), dimension(:), pointer :: surfs  => null()
      real(rv), dimension(:,:,:), pointer :: values => null()
      logical :: Verbose
      ! Executable code
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_begin ( 'FillUtils_1.NamedQtyFromFile', key )
      verbose = ( switchDetail(switches, 'fill') > -1 )
      call GetHDF5DSDims ( MLSFile, name, DIMS )
      dimInts = max(dims, int(1,hsize_t))
      select case (lowercase(fileType))
      case ('l2aux')
        if ( spread ) call Announce_Error ( key, no_Error_Code, &
          &   'Unable to use spread when filling from l2aux file' )
        call ReadL2AUXData ( MLSFile, name, quantity%template%quantityType, l2aux )
        call FromL2Aux( quantity, l2aux, status )
        call DestroyL2AUXContents ( l2aux )
      case ('swath', 'l2gp')
        if ( spread ) call Announce_Error ( key, no_Error_Code, &
          &   'Unable to use spread when filling from l2gp file' )
        call ReadL2GPData( MLSFile, name, l2gp )
        call FromL2GP ( quantity, l2gp, .false., -1, status, .true., .false. )
        call DestroyL2GPContents ( l2gp )
      case default ! E.g., 'hdf'
        call Allocate_test ( values, dimInts(1), dimInts(2), dimInts(3), 'values read from file', ModuleName )
        call loadFromHDF5DS ( MLSFile, &
          & trim(Name), values ) ! quantity%values )
        if ( verbose ) then
          call outputNamedValue( 'name', trim(name) )
          call outputNamedValue( 'spread', spread )
          call outputNamedvalue( 'shape(quantity%values)', shape(quantity%values) )
          call outputNamedvalue( 'shape(values)', shape(values) )
        endif
        if ( spread ) then
          quantity%values = values(1,1,1)
        elseif ( interpolate .and. associated(quantity%template%surfs) ) then
          ! Must we interpolate according to surface heights?
          call h5dOpen_f ( MLSFile%fileID%f_id, trim(name), &
            & MLSFile%fileID%sd_id, status)
          if ( status /= 0 ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & 'Unable to open sd to read heights attribute in hdf file', &
                  & MLSFile=MLSFile )
          end if
          call GetHDF5Attribute( MLSFile, 'noSurfs', noSurfs )
          call Allocate_test ( surfs, noSurfs, 'surfs read from file', ModuleName )
          call GetHDF5Attribute( MLSFile, 'surfs', surfs )
          if ( all(surfs == quantity%template%surfs(:,1)) ) then
            ! No need, the heights are all the same
            call fillMyInstances( quantity%values, values )
          elseif ( noSurfs == dimInts(1)*dimInts(2) ) then
            call FromProfile( quantity, surfs, &
              & real( &
              & reshape( values(:,:,1), (/ dimInts(1)*dimInts(2) /) ), &
              & rv), &
              & instancesNode, globalUnit, ptan )
          else
            call Announce_Error ( key, no_Error_Code, &
          &   'Unable to interpolate using this file data set ; noSurfs ' // &
          &   'does not match dims(1)*dims(2)' )
          endif
          call DeAllocate_test ( surfs, 'surfs read from file', ModuleName )
          call h5dClose_f ( MLSFile%fileID%sd_id, status)
          if ( status /= 0 ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & 'Unable to close sd to read heights attribute in hdf file', &
                  & MLSFile=MLSFile )
          end if
          MLSFile%fileID%sd_id = 0
        elseif ( size(quantity%values(:,1)) == size(values(:,:,1) ) ) then
          do instance = 1, size(quantity%values(1,:))
            quantity%values(:,instance) = &
              & reshape( values(:,:,1), (/ dimInts(1)*dimInts(2) /) )
          enddo
        elseif ( size(quantity%values) /= size(values ) ) then
          call Announce_Error ( 0, no_Error_Code, &
            'sizes of read and filled datasets dont match' )
        else
          quantity%values = reshape( values, (/ dimInts(1)*dimInts(2), dimInts(3) /) )
        endif
        call DeAllocate_test ( values, 'values read from file', ModuleName )
      end select
      if ( toggle(gen) .and. levels(gen) > 2 ) &
        & call trace_end ( 'FillUtils_1.NamedQtyFromFile' )
    contains
      subroutine FillMyInstances( values1, values2 )
        ! Handle the following cases:
        ! (1) values2 has only 1 instance
        ! (2) values2 has the same number of instances as values1
        ! (3) values2 has a different number of instances from values1
        ! Args:
        real(rv), dimension(:,:), pointer :: VALUES1
        real(rv), dimension(:,:,:), pointer :: values2
        ! Local variables
        integer :: instance
        integer :: noInst1  ! Number of instances in values1
        integer :: noInst2  ! Number of instances in values2
        integer, dimension(3) :: shp2
        ! Executable code
        noInst1 = size(values1, 2)
        noInst2 = size(values2, 3)
        shp2 = shape(values2)
        if ( noInst2 == 1 .or. homogeneous ) then
          do instance = 1, noInst1
            values1(:,instance) = &
              & reshape( values2(:,:,1), (/ shp2(1)*shp2(2) /) )
          enddo
        elseif( noInst1 == noInst2 ) then
          do instance = 1, noInst1
            values1(:,instance) = &
              & reshape( values2(:,:,instance), (/ shp2(1)*shp2(2) /) )
          enddo
        else
          do instance = 1, noInst1
            values1(:,instance) = &
              & reshape( values2(:,:,1), (/ shp2(1)*shp2(2) /) )
          enddo
        endif
      end subroutine FillMyInstances
    end subroutine NamedQtyFromFile
    
  ! --------- WhichSurfaceIsHeight ------
  function WhichSurfaceIsHeight ( node, Heights ) result ( surface )
    ! Args
    integer, intent(in)                  :: node  ! Tree node
    real(r8), dimension(:), intent(in)   :: heights ! array of heights
    integer                              :: surface
    ! Internal variables
    real(r8) :: HEIGHT                ! The height to consider
    integer, dimension(2) :: unitAsArray ! Unit for value given
    logical :: unitsError
    real (r8), dimension(2) :: valueAsArray ! Value given
    ! Executable
    if ( node > 0 ) then
      if ( nsons ( node ) /= 2 ) call Announce_Error ( node, no_error_code, &
        & 'Only one height can be supplied for explicit fill' )
      call expr_check ( subtree(2,node) , unitAsArray, valueAsArray, &
        & (/PHYQ_Pressure/), unitsError )
      if ( unitsError ) call Announce_error ( node, wrongUnits, &
        & extraInfo=(/unitAsArray(1), PHYQ_Pressure/) )
      height = - log10 ( valueAsArray(1) )
      call Hunt ( heights, height, surface, nearest=.true. )
    else
      surface = 0
    endif
  end function WhichSurfaceIsHeight

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module FillUtils_1
!=============================================================================

!
! $Log$
! Revision 2.75  2013/05/21 01:47:34  vsnyder
! Don't look at QUITNOW until you're sure it's there
!
! Revision 2.74  2013/05/17 00:52:01  pwagner
! May constrain Transfer command command to quantitynames by r/t Boolean; r/t sep now achar(0)
!
! Revision 2.73  2013/04/05 23:20:47  pwagner
! Requires verbose settings for extra printing
!
! Revision 2.72  2013/02/01 23:42:25  vsnyder
! Use CreateVectorValue instead of Allocate_Test
!
! Revision 2.71  2013/01/14 21:22:51  pwagner
! Changed chiSqRatio when skipped to 999
!
! Revision 2.70  2013/01/02 21:41:31  pwagner
! Added derivative method to Fill command; Transfer can do Fill methods, too
!
! Revision 2.69  2012/12/20 01:05:26  vsnyder
! Provide coherent-to-minor frame interpolation
!
! Revision 2.68  2012/11/14 00:59:23  pwagner
! Use dimList for choosing which of {csi} to average over
!
! Revision 2.67  2012/11/08 23:21:51  pwagner
! dimList field lets us specifiy whether to shift by [c,s,i] during manipulate
!
! Revision 2.66  2012/11/05 19:04:28  pwagner
! Most Fill methods need not handle dontMask themselves
!
! Revision 2.65  2012/10/22 18:15:47  pwagner
! Many Subset operations now available in Fill
!
! Revision 2.64  2012/10/09 00:48:31  pwagner
! New ignoreTemplate, changed force meaning in Fill
!
! Revision 2.63  2012/08/16 17:58:00  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.62  2012/08/08 19:57:06  vsnyder
! Use a matrix block instead of a real array in FillCovariance
!
! Revision 2.61  2012/07/31 00:49:02  vsnyder
! Use DestroyVectoryQuantityMask abstraction
!
! Revision 2.60  2012/06/07 22:46:03  pwagner
! Turn off warnings about non-monotonic Ptan in FromProfile
!
! Revision 2.59  2012/06/07 00:53:25  vsnyder
! More, and more consistent, -g tracing
!
! Revision 2.58  2012/05/08 17:51:58  pwagner
! Disabled another debugging bloc
!
! Revision 2.57  2012/04/26 23:30:04  pwagner
! Monotonize ptan heights when needed; (should we sort instead?)
!
! Revision 2.56  2012/04/20 01:06:03  pwagner
! Changes to bin, ChiSqRatio, and Interpolating to remove bugs, hopefully adding none
!
! Revision 2.55  2012/03/12 17:11:46  pwagner
! 'ln' same as 'log' in manipulate fills; new 'cgrid' switch
!
! Revision 2.54  2012/02/24 21:22:02  pwagner
! DirectRead may /interpolate vertically
!
! Revision 2.53  2012/02/16 22:39:02  pwagner
! Dont automatically in AutoFillVector
!
! Revision 2.52  2012/01/25 01:18:05  pwagner
! Fixed most bugs in QtyFromFile and VectorFromFile; now can select by attributes
!
! Revision 2.51  2012/01/18 20:38:36  vsnyder
! Add 'matrix and vector inconsistent' error message
!
! Revision 2.50  2011/12/15 01:50:26  pwagner
! Added sdName and /spread fields to DirectRead
!
! Revision 2.49  2011/11/30 21:30:00  pwagner
! Fixed bug affecting skipped retrievals
!
! Revision 2.48  2011/11/04 23:46:45  pwagner
! Fixed bug added with last commit
!
! Revision 2.47  2011/11/04 00:24:39  pwagner
! Added procedures to auto special qties in a vector, fill a qty from a file, and fill a vector from a file
!
! Revision 2.46  2011/07/20 00:53:40  pwagner
! Fixed bug added in rev2.43
!
! Revision 2.45  2011/06/16 20:51:24  vsnyder
! Make Announce_Error codes public
!
! Revision 2.44  2011/04/20 16:46:00  pwagner
! Fixed solecism that NAG complained about
!
! Revision 2.43  2011/04/13 00:28:15  pwagner
! manipulation can now accept digits
!
! Revision 2.42  2011/04/07 23:35:31  pwagner
! Fixed bug in handling manipulation='(c - c*a)*b'
!
! Revision 2.41  2011/03/08 18:29:22  pwagner
! Added shift,slip,chaannel,surface,instance functions to manipulate fills
!
! Revision 2.40  2010/07/06 16:05:37  pwagner
! Better error checking in MaanipulateVectors
!
! Revision 2.39  2010/07/01 00:49:33  pwagner
! Transfer between vectors may now also manipulate
!
! Revision 2.38  2010/06/18 16:48:34  pwagner
! Corrected error that prevented filling radiances qty from l1b
!
! Revision 2.37  2010/05/24 16:33:52  honghanh
! Merge changes in version 2.35 and 2.36
!
! Revision 2.35  2010/05/19 23:06:45  pwagner
! Shorten most Fill routine names
!
! Revision 2.34  2010/04/28 16:23:52  pwagner
! May specify instances range in explicit Fill
!
! Revision 2.33  2010/04/22 23:36:00  pwagner
! May fill num rads/MIF as a percentage
!
! Revision 2.32  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.31  2009/12/14 18:35:51  pwagner
! Dont crash in FromGrid if Grid is empty
!
! Revision 2.30  2009/10/27 22:14:24  pwagner
! Compiles with new api for Dump vector quantity
!
! Revision 2.29  2009/09/01 17:14:02  pwagner
! Reduce severity of profile mismatch in FromL2GP to permit filling 1d sids with truth
!
! Revision 2.28  2009/08/24 20:14:11  pwagner
! May Fill H2O precision from RHI precision
!
! Revision 2.27  2009/07/21 20:34:56  pwagner
! chi^2 ratio nay hold values for iterations prior to final
!
! Revision 2.26  2009/07/10 20:56:52  pwagner
! Fixed bug affecting manipulations like 'a+b+c'
!
! Revision 2.25  2009/06/30 15:19:05  pwagner
! Fixed bug in Explicit
!
! Revision 2.24  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.23  2009/05/13 20:41:55  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.22  2009/05/01 23:44:40  pwagner
! Restored nearly all except for conversions between RHi and H2O
!
! Revision 2.21  2009/04/30 22:13:29  pwagner
! name of bit in MaskVectorQty and isVectorQtyMasked now mandatory
!
! Revision 2.20  2009/04/30 20:15:01  pwagner
! Another fix to masking bit miscues in RHI..
!
! Revision 2.19  2009/04/29 23:11:04  pwagner
! Explicit can Fill from an optional extraQuantity
!
! Revision 2.18  2009/04/28 20:03:49  pwagner
! Consider only mask of quantity being filled, not sources in RHi
!
! Revision 2.17  2009/04/16 21:55:23  pwagner
! /exact keyword in status Fill to fix radiance bug
!
! Revision 2.16  2009/04/13 20:45:57  pwagner
! heightRange in explicit Fill can fill above or below specified height
!
! Revision 2.15  2009/03/05 18:38:32  pwagner
! May specifiy height, channel with explicit Fill
!
! Revision 2.14  2008/10/15 16:37:28  pwagner
! Let precisions explicitly set negative also mask radiances
!
! Revision 2.13  2008/09/24 16:46:09  livesey
! Changed ptan from optional to pointer in fill from profile
!
! Revision 2.12  2008/09/20 00:03:00  pwagner
! Added print statement to not_used_here
!
! Revision 2.11  2008/09/16 22:30:19  pwagner
! Use optional arg ptan when source is profile, vector
!
! Revision 2.10  2008/08/14 20:58:40  pwagner
! /interpolate now possible field in Transfer command
!
! Revision 2.9  2008/08/06 17:27:47  pwagner
! Fill by manipulation now respects mask better
!
! Revision 2.8  2008/06/06 21:02:49  michael
! added fill method uncompressRadiance
!
! Revision 2.7  2008/04/26 00:39:56  livesey
! Added total power stuff
!
! Revision 2.6  2008/04/11 01:17:09  livesey
! Added uncompressRadiance fill
!
! Revision 2.5  2008/01/07 21:43:03  pwagner
! Fixed bug regarding unrecognized ops '<' and '>' in catTwoOperands
!
! Revision 2.4  2007/11/05 18:41:20  pwagner
! expr with 'c' unambiguous (I hope); '^' power op added to expr
!
! Revision 2.3  2007/11/01 23:33:59  pwagner
! rewrite to permit functions and algebra in same manipulation; needs more testing
!
! Revision 2.2  2007/10/04 20:43:12  vsnyder
! Remove unused symbols
!
! Revision 2.1  2007/09/27 21:59:00  pwagner
! First commit
!
