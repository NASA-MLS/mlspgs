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
module Output_UARS_L1B
!=============================================================================

  use MLS_DataProducts, only: DataProducts_T, Deallocate_DataProducts
  use MLSAuxData, only: Build_MLSAuxData, CreateGroup_MLSAuxData
  use MLSSTRINGS, only: CompressString

  implicit NONE

  private

  integer, parameter :: Int_Fill = -1

  public :: OutputL1B_OA, OutputL1B_Rad

  logical, parameter, private :: Deebug = .false.

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  subroutine OutputL1B_rad ( sdId, noMAF, counterMAF, rad, prec, gird_time, &
                           & Limb_Hdr )

    use MLSKinds, only: R8
    use Rad_File_Contents, only: Limb_Hdr_t

   ! This subroutine writes a MAF's worth of data to the L1BRad File

    integer, intent(in) :: sdId, noMAF, counterMAF
    real, intent(in) :: rad(15,6,32), prec(15,6,32)
    real(r8), intent(in) :: gird_time
    type(limb_hdr_t), intent(in) :: Limb_Hdr

    TYPE (DataProducts_T) :: dataset
    integer, dimension(3) :: dims = (/15, 32, 1/)
    integer :: bank, status, bits, spindx

    real, parameter :: Rad_Fill = -999.99, Rad_Err_Fill = -1.0

    character(len=40) :: name
    character(len=4), parameter :: species(6) = (/ &
         "PT  ", "ClO " , "H2O2", "O3  ", "H2O ", "O3  " /)
    character(len=3), parameter :: band(6) = (/ &
         'B1F', 'B2F', 'B3F', 'B4F', 'B5F', 'B6F' /)
    character(len=2), parameter :: radiometer(3) = (/ 'R1', 'R2', 'R3' /)
    character(len=3), parameter :: freq(6) = (/ &
         '63 ', '205', '205', '205', '183', '183' /)
    character(len=2), parameter :: switch(0:2) = (/ 'S0', 'S1', 'S2' /)
    integer :: mask07, bank3, bank6, swpos
    data mask07 / z'07'/  ! lower 3 bits for Filter Bank 3

    call Deallocate_DataProducts (dataset)

    allocate (dataset%Dimensions(1), stat=status)

    dataset%name      = 'counterMAF'
    dataset%data_type = 'integer'
    dataset%Dimensions(1) = 'MAF'
    call Build_MLSAuxData (sdId, dataset, counterMAF, fill_value=INT_FILL, &
         lastIndex=noMAF)

    dataset%name      = 'MAFStartTimeGIRD'
    dataset%data_type = 'double'
    call Build_MLSAuxData (sdId, dataset, gird_time, fill_value=-1.0d0, &
         lastIndex=noMAF)

    deallocate (dataset%Dimensions, stat=status)

    allocate (dataset%Dimensions(3), stat=status)
    dataset%data_type = 'real'
    dataset%Dimensions(1) = 'chanFB'
    dataset%Dimensions(2) = 'MIF'
    dataset%Dimensions(3) = 'MAF'

    bits = ichar(limb_hdr%band_bank)

    bank3 = iand (bits, mask07)  ! what's hooked up to filter bank 3
    if (btest (bits, 3)) then    ! what's hooked up to filter bank 6
       bank6 = 3
    else
       bank6 = 6
    end if

    do bank = 1, 6
       swpos = 0
       if (bank == 3) then
          spindx = bank3
          if (bank3 /= 3) swpos = 1
       else if (bank == 6) then
          spindx = bank6
          if (bank6 /= 3) swpos = 2
       else
          spindx = bank
       end if
       if (bank == 1) then
          name = radiometer(1)
       else if (bank >= 2 .and. bank <= 4) then
          name = radiometer(2)
       else
          name = radiometer(3)
       end if

       name = name(1:len_trim(name)) // ':' // freq(bank) // '.' // &
            band(bank) // ':' // species(spindx) // '.' // switch(swpos)
       name = name(1:len_trim(name)) // '.FB15-' // char(bank+48)
       name = CompressString (name) 

       dataset%name = TRIM(name)
       call Build_MLSAuxData (sdId, dataset, rad(:,bank,:), &
            lastIndex=noMAF, dims=dims, fill_value=RAD_FILL)
       dataset%name = TRIM(name)//' precision'
       call Build_MLSAuxData (sdId, dataset, prec(:,bank,:), &
            lastIndex=noMAF, dims=dims, fill_value=RAD_ERR_FILL)
    end do

  end subroutine OutputL1B_rad

  subroutine OutputL1B_OA (sdId, noMAF, counterMAF, Limb_OA, EMLS_OA )

    ! This subroutine writes a MAF's worth of data to the L1B OA file

    use HDF5, only:  H5GClose_f, H5GOpen_f
    use MLSHDF5, only: MakeHDF5Attribute
    use MLSKinds, only: R8
    use oa_file_contents, only: emls_oa_t
    use rad_file_contents, only: limb_oa_t
 
    ! Arguments:
    integer, intent(in) :: sdId, noMAF, counterMAF
    type(limb_oa_t), intent(in) :: Limb_OA
    type(emls_oa_t), intent(in) :: EMLS_OA

    integer :: Grp_Id
    integer :: status

    type (DataProducts_T) :: dataset

!------------------------------------------------------------------------------
    call Deallocate_DataProducts (dataset)

    if ( deebug ) print *, 'MAF: ', noMAF, limb_oa%ref_lat, limb_oa%ref_time

! 2-d first:

    allocate (dataset%Dimensions(2), stat=status)

    dataset%data_type = 'real'
    dataset%Dimensions(1) = 'MIF'
    dataset%Dimensions(2) = 'MAF'

    dataset%name = 'sc/OrbIncl'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_orbincl, lastIndex=noMAF)
    dataset%name = 'sc/Lon'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_lon, lastIndex=noMAF)
    dataset%name = 'sc/GeocLat'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_geoclat, lastIndex=noMAF)
    dataset%name = 'sc/GeodLat'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_geodlat, lastIndex=noMAF)
    dataset%name = 'GHz/GeocLat'
    call Build_MLSAuxData (sdId, dataset, emls_oa%geoc_lat, lastIndex=noMAF)
    dataset%name = 'GHz/GeodLat'
    call Build_MLSAuxData (sdId, dataset, emls_oa%geod_lat, lastIndex=noMAF)
    dataset%name = 'GHz/GeodAngle'
    call Build_MLSAuxData (sdId, dataset, emls_oa%geod_ang, lastIndex=noMAF)
    dataset%name = 'GHz/Lon'
    call Build_MLSAuxData (sdId, dataset, emls_oa%lon, lastIndex=noMAF)
    dataset%name = 'GHz/SolarTime'
    call Build_MLSAuxData (sdId, dataset, emls_oa%solartime, lastIndex=noMAF)
    dataset%name = 'GHz/SolarZenith'
    call Build_MLSAuxData (sdId, dataset, emls_oa%solarzenith, lastIndex=noMAF)
    dataset%name = 'GHz/azimAngle'
    call Build_MLSAuxData (sdId, dataset, emls_oa%azimAngle, lastIndex=noMAF)
    dataset%name = 'GHz/scanAngle'
    call Build_MLSAuxData (sdId, dataset, emls_oa%scanAngle, lastIndex=noMAF)

    dataset%data_type = 'double'
    dataset%name = 'GHz/GeodAlt'
    call Build_MLSAuxData (sdId, dataset, emls_oa%geod_alt, lastIndex=noMAF)
    dataset%name = 'GHz/GeocAlt'
    call Build_MLSAuxData (sdId, dataset, emls_oa%geoc_alt, lastIndex=noMAF)
    dataset%name = 'GHz/LosAngle'
    call Build_MLSAuxData (sdId, dataset, emls_oa%LosAngle, lastIndex=noMAF)
    dataset%name = 'GHz/LosVel'
    call Build_MLSAuxData (sdId, dataset, emls_oa%LosVel, lastIndex=noMAF)
    dataset%name = 'GHz/OrbY'
    call Build_MLSAuxData (sdId, dataset, emls_oa%OrbY, lastIndex=noMAF)
    dataset%name = 'GHz/ECRtoFOV'
    call Build_MLSAuxData (sdId, dataset, emls_oa%ECRtoFOV, lastIndex=noMAF)
    dataset%name = 'sc/GeocAlt'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_geocalt, lastIndex=noMAF)
    dataset%name = 'sc/GeodAlt'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_geodalt, lastIndex=noMAF)

    dataset%data_type = 'integer'
    dataset%name = 'GHz/BO_stat'
    call Build_MLSAuxData (sdId, dataset, emls_oa%bo_stat, lastIndex=noMAF)

    allocate (dataset%Dimensions(1), stat=status)
    dataset%data_type = 'double'
    dataset%Dimensions(1) = 'MAF'
    dataset%name = 'MAFStartTimeTAI'
    call Build_MLSAuxData (sdId, dataset, emls_oa%MAFStartTimeTAI, &
         lastIndex=noMAF)
    dataset%name = 'sc/MIF_TAI'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_MIF_TAI, lastIndex=noMAF)

    dataset%data_type = 'integer'
    dataset%name      = 'counterMAF'
    call Build_MLSAuxData (sdId, dataset, counterMAF, fill_value=INT_FILL, &
         lastIndex=noMAF)

    call Deallocate_DataProducts (dataset)

! 3-d next:

    deallocate (dataset%Dimensions, stat=status)
    allocate (dataset%Dimensions(3), stat=status)

    dataset%data_type = 'double'
    dataset%Dimensions(1) = 'xyz'
    dataset%Dimensions(2) = 'MIF'
    dataset%Dimensions(3) = 'MAF'
    dataset%name      = 'GHz/ECI'
    call Build_MLSAuxData (sdId, dataset, emls_oa%ECI, lastIndex=noMAF)
    dataset%name      = 'GHz/ECR'
    call Build_MLSAuxData (sdId, dataset, emls_oa%ECR, lastIndex=noMAF)
    dataset%name      = 'sc/ECI'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_ECI, lastIndex=noMAF)
    dataset%name      = 'sc/ECR'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_ECR, lastIndex=noMAF)
    dataset%name      = 'sc/VelECI'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_VelECI, lastIndex=noMAF)
    dataset%name      = 'sc/VelECR'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_VelECR, lastIndex=noMAF)
    dataset%name      = 'sc/ypr'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_ypr, lastIndex=noMAF)
    dataset%name      = 'sc/yprRate'
    call Build_MLSAuxData (sdId, dataset, emls_oa%sc_ypr_rate, lastIndex=noMAF)

    call Deallocate_DataProducts (dataset)

    ! Now attributes for some component units
    call h5gopen_f ( sdId, '/GHz', grp_id, status )
    call MakeHDF5Attribute ( grp_id, 'GeocLat_units',       'deg'  )
    call MakeHDF5Attribute ( grp_id, 'GeodLat_units',       'deg'  )
    call MakeHDF5Attribute ( grp_id, 'GeodAng_units',       'deg'  )
    call MakeHDF5Attribute ( grp_id, 'Lon_units',           'deg'  )
    call MakeHDF5Attribute ( grp_id, 'GeocAlt_units',       'm'    )
    call MakeHDF5Attribute ( grp_id, 'GeodAlt_units',       'm'    )
    call MakeHDF5Attribute ( grp_id, 'LosAngle_units',      'deg'  )
    call MakeHDF5Attribute ( grp_id, 'LosVel_units',        'm/s'  )
    call MakeHDF5Attribute ( grp_id, 'OrbY_units',          'm'    )
    call MakeHDF5Attribute ( grp_id, 'ECI_units',           'm'    )
    call MakeHDF5Attribute ( grp_id, 'ECR_units',           'm'    )
    call h5gclose_f ( grp_id, status )
    call h5gopen_f ( sdId, '/sc', grp_id, status )
    call MakeHDF5Attribute ( grp_id, 'OrbIncl_units',       'deg'  )
    call MakeHDF5Attribute ( grp_id, 'Sc_Lon_units',        'deg'  )
    call MakeHDF5Attribute ( grp_id, 'Sc_GeocLat_units',    'deg'  )
    call MakeHDF5Attribute ( grp_id, 'Sc_GeodLat_units',    'deg'  )
    call MakeHDF5Attribute ( grp_id, 'Sc_GeocAlt_units',    'm'    )
    call MakeHDF5Attribute ( grp_id, 'Sc_GeodAlt_units',    'm'    )
    call MakeHDF5Attribute ( grp_id, 'ECI_units',           'm'    )
    call MakeHDF5Attribute ( grp_id, 'ECR_units',           'm'    )
    call MakeHDF5Attribute ( grp_id, 'VelECI_units',        'm/s'  )
    call MakeHDF5Attribute ( grp_id, 'VelECR_units',        'm/s'  )
    call h5gclose_f ( grp_id, status )

  end subroutine OutputL1B_OA

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Output_UARS_L1B

! $Log$
! Revision 1.2  2014/12/11 00:48:51  vsnyder
! Move external procedures into modules.  Add copyright and CVS lines.
! Compute MIF geolocation (except height) for SC.  Compute MIF-resolved
! SC velocity.  Some cannonball polishing.
!
