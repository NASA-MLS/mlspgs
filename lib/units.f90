! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module UNITS
! Provide several units and unit conversion constants.
! Initialize the declaration table with the unit names and scales.
  use DECLARATION_TABLE, only: DECLARE, UNITS_NAME
  use INTRINSIC ! "units" type literals, beginning with L_,
                ! Abstract physical quantities beginning with PHYQ_
  use MLSCommon, only: R8
  use TREE, only: NULL_TREE

  implicit NONE
  public

  real(r8), parameter :: Ln10 = 2.302585092994045684017991454684364207601
  real(r8), parameter :: Pi = 3.141592653589793238462643383279502884197
  real(r8), parameter :: Deg2Rad = Pi/180.0_r8 ! Degrees-to-Radians
  real(r8), parameter :: Rad2Deg = 180.0_r8/Pi ! Radians-to-Degrees
  real(r8), parameter :: Omega = 7.292115D-5 ! Angular velocity of earth

  private :: DECLARE, NULL_TREE, UNITS_NAME ! Get them from the source
  ! But it's OK to get PHYQ_... from here

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public procedures     =============================
! ---------------------------------------------------  INIT_UNITS  -----
  subroutine INIT_UNITS
    ! Put scale factors for units into the declaration table.

    ! If the scale factor is negative, it is subtracted instead of
    ! multiplied.  So far, this is used only to convert Celsius to
    ! Kelvin

    ! Call this procedure AFTER "init_tables_module%init_tables"
    ! We can't call it from there because there would be a circular
    ! dependence.  We can't put it there because a circular dependence
    ! between "init_tables_module" and "declaration_table" would
    ! result.

    call declare_unit ( l_dimensionless, 1.0d0, phyq_dimensionless )
    call declare_unit ( l_dimless, 1.0d0, phyq_dimensionless )
    call declare_unit ( l_dl, 1.0d0, phyq_dimensionless )

    call declare_unit ( l_m, 1.0d0, phyq_length )
    call declare_unit ( l_meters, 1.0d0, phyq_length )
    call declare_unit ( l_km, 1000.0d0, phyq_length )

    call declare_unit ( l_s, 1.0d0, phyq_time )
    call declare_unit ( l_seconds, 1.0d0, phyq_time )
    call declare_unit ( l_minutes, 60.0d0, phyq_time )
    call declare_unit ( l_hours, 60.0d0 * 60.0d0, phyq_time )
    call declare_unit ( l_days, 60.0d0 * 60.0d0 * 24.0d0, phyq_time )

    call declare_unit ( l_mb, 1.0d0, phyq_pressure )
    call declare_unit ( l_hpa, 1.0d0, phyq_pressure )
    call declare_unit ( l_pa, 1.0d-2, phyq_pressure )

    call declare_unit ( l_k, 1.0d0, phyq_temperature )
    call declare_unit ( l_c, -273.15d0, phyq_temperature )

    call declare_unit ( l_vmr, 1.0d0, phyq_vmr )
    call declare_unit ( l_ppmv, 1.0d-6, phyq_vmr )
    call declare_unit ( l_ppbv, 1.0d-9, phyq_vmr )
    call declare_unit ( l_pptv, 1.0d-12, phyq_vmr )

    call declare_unit ( l_deg, 1.0d0, phyq_angle )
    call declare_unit ( l_degrees, 1.0d0, phyq_angle )
    call declare_unit ( l_radians, 45.0d0 / atan(1.0d0), phyq_angle )
    call declare_unit ( l_rad, 45.0d0 / atan(1.0d0), phyq_angle )
    call declare_unit ( l_orbits, 360.0d0, phyq_angle )

    call declare_unit ( l_maf, 1.0d0, phyq_mafs )
    call declare_unit ( l_mafs, 1.0d0, phyq_mafs )

    call declare_unit ( l_mif, 1.0d0, phyq_mifs )
    call declare_unit ( l_mifs, 1.0d0, phyq_mifs )

    call declare_unit ( l_mhz, 1.0d0, phyq_frequency )
    call declare_unit ( l_ghz, 1.0d3, phyq_frequency )
    call declare_unit ( l_thz, 1.0d6, phyq_frequency )
    call declare_unit ( l_khz, 1.0d-3, phyq_frequency )
    call declare_unit ( l_hz, 1.0d-6, phyq_frequency )

    call declare_unit ( l_zeta, 1.0d0, phyq_zeta )
    call declare_unit ( l_logp, 1.0d0, phyq_zeta )
!   this is = 1/meters 
    call declare_unit ( l_invm, 1.0d-3, phyq_extinction )
    call declare_unit ( l_invkm, 1.0d0, phyq_extinction )
!    call declare_unit ( l_extinction, 1.0d0, phyq_extinction )
!   this is = 1Pa*1000*1sec^2/1meter^4
    call declare_unit ( l_icedensity, 1.0d1, phyq_icedensity )
!   this is 1DU ( = 2.687e20 molecules/m^2)
    call declare_unit ( l_DobsonUnits, 1.0d0, phyq_DobsonUnits )

  contains
    subroutine DECLARE_UNIT ( NAME, VALUE, PHYS_UNIT )
      integer, intent(in) :: NAME
      integer, intent(in) :: PHYS_UNIT
      double precision, intent(in) :: VALUE
    ! Declare "name" to be a unit of abstract kind "phys_unit" (e.g. length)
    ! and a scale to that kind of "value"
      call declare ( lit_indices(name), value, units_name, phys_unit, &
                     null_tree )
    end subroutine DECLARE_UNIT
  end subroutine INIT_UNITS
end module UNITS

! $Log$
! Revision 2.12  2001/12/06 23:45:34  livesey
! Moved Omega into here. Might be time to have a physical constants
! module somehere?
!
! Revision 2.11  2001/11/08 00:12:31  livesey
! Commented out extinction as it's now a molecule, and so not known yet.
! Perhaps we'll sort this out later.
!
! Revision 2.10  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.9  2001/07/10 23:47:01  jonathan
! added phyq_icedensity, paul/jonathan
!
! Revision 2.8  2001/04/26 02:35:38  vsnyder
! Fix up CVS stuff
!
! Revision 2.7  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.6  2001/04/09 20:59:35  vsnyder
! Add C (for Celsius) unit and l_c name for it
!
! Revision 2.5  2001/03/27 18:42:08  vsnyder
! Insert mathematical constants
!
! Revision 2.4  2001/03/21 23:28:43  livesey
! Fixed bug in l_zeta and l_logp, had inappropriate scaling factors
!
! Revision 2.3  2000/10/27 21:55:12  pwagner
! Never should have left lib
!
! Revision 1.1  2000/10/27 20:30:35  pwagner
! moved to test_cf_parser
!
! Revision 2.1  2000/10/11 18:23:40  vsnyder
! Insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:52  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
