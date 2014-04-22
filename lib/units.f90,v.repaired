! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module UNITS

! Provide several units and unit conversion constants.
! Initialize the declaration table with the unit names and scales.

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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

    use Constants, only: Rad2Deg
    use Intrinsic ! "units" type literals, beginning with L_,
                  ! Abstract physical quantities beginning with PHYQ_
    use Toggles, only: Con, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    integer :: Me = -1 ! String index for tracing

    call trace_begin ( me, 'INIT_UNITS', cond=toggle(con) .and. levels(con) > 1 )

    call declare_unit ( l_dl, 1.0d0, phyq_dimensionless )
    call declare_unit ( l_dimless, 1.0d0, phyq_dimensionless )
    call declare_unit ( l_dimensionless, 1.0d0, phyq_dimensionless )

    call declare_unit ( l_m, 1.0d0, phyq_length )
    call declare_unit ( l_meters, 1.0d0, phyq_length )
    call declare_unit ( l_km, 1000.0d0, phyq_length )

    call declare_unit ( l_us, 1.0d-6, phyq_time )
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
    call declare_unit ( l_radians, rad2deg, phyq_angle )
    call declare_unit ( l_rad, rad2deg, phyq_angle )
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

    call declare_unit ( l_ms, 1.0d0, phyq_velocity )
    call declare_unit ( l_kms, 1.0d3, phyq_velocity )

    call declare_unit ( l_zeta, 1.0d0, phyq_zeta )
    call declare_unit ( l_logp, 1.0d0, phyq_zeta )

!   this is = 1/meters 
    call declare_unit ( l_invm, 1.0d-3, phyq_extinction )
    call declare_unit ( l_invkm, 1.0d0, phyq_extinction )

!    call declare_unit ( l_extinction, 1.0d0, phyq_extinction )
!    call declare_unit ( l_extinctionV2, 1.0d0, phyq_extinction )
!    call declare_unit ( l_MIFextinction, 1.0d0, phyq_extinction )
!    call declare_unit ( l_MIFextinctionV2, 1.0d0, phyq_extinction )

!   this is = 1Pa/1000*1sec^2/1meter^2 = 1gm/meter^3
    call declare_unit ( l_icedensity, 1.0d0, phyq_icedensity )
    call declare_unit ( l_iwc, 1.0d0, phyq_icedensity )
    call declare_unit ( l_gm3, 1.0d0, phyq_icedensity )
    call declare_unit ( l_mgm3, 1.0d-3, phyq_icedensity )

!   this is 1DU ( = 2.687e20 molecules/m^2)
!   (but we will use molecules/cm^2 as the default)
    call declare_unit ( l_DobsonUnits, 2.687d16, phyq_colmabundance )
    call declare_unit ( l_DU, 2.687d16, phyq_colmabundance )
    call declare_unit ( l_molcm2, 1.0d0, phyq_colmabundance )

    call declare_unit ( l_pctrhi, 1.0d0, phyq_pctrhi )

    call declare_unit ( l_gauss, 1.0d0, phyq_gauss )

    call declare_unit ( l_profiles, 1.0d0, phyq_profiles )

    call trace_end ( 'INIT_UNITS', cond=toggle(con) .and. levels(con) > 1 )

  contains

    subroutine DECLARE_UNIT ( NAME, VALUE, PHYS_UNIT )
      use DECLARATION_TABLE, only: DECLARE, PHYS_UNIT_NAME, UNITS_NAME
      use Intrinsic, only: PHYQ_Indices
      integer, intent(in) :: NAME
      integer, intent(in) :: PHYS_UNIT
      double precision, intent(in) :: VALUE
      integer :: Me = -1 ! String index for tracing

      call trace_begin ( me, 'DECLARE_UNIT', cond=toggle(con) .and. levels(con) > 8 )
    ! Declare "name" to be a unit of abstract kind "phys_unit" (e.g. length)
    ! and a scale to that kind of "value"
      !              string             value  type        units
      call declare ( lit_indices(name), value, units_name, phys_unit, &
      !              tree
                     name )
      !              string                   value  type
      call declare ( phyq_indices(phys_unit), value, phys_unit_name, &
      !              units      tree
                     phys_unit, name )
      call trace_end ( 'DECLARE_UNIT', stringIndex=lit_indices(name), &
        & cond=toggle(con) .and. levels(con) > 8 )
    end subroutine DECLARE_UNIT

  end subroutine INIT_UNITS

! ----------------------------------------------------  Base_Unit  -----
  integer function Base_Unit ( PHYS_Unit )
    ! Get the base unit corresponding to a physical unit.
    ! The base unit is one with scale == 1.0.
    ! The result is the base unit's lit index, not it's string index.
    use Declaration_Table, only: Decls, Get_Decl, Phys_Unit_Name
    use Intrinsic, only: L_Dimensionless, PHYQ_Indices
    integer, intent(in) :: PHYS_Unit
    type(decls) :: Decl
    base_unit = l_dimensionless ! for want of a better default
    decl = get_decl ( phyq_indices(phys_unit), type=phys_unit_name, value=1.0d0 )
    if ( decl%type == phys_unit_name ) base_unit = decl%tree
  end function Base_Unit

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module UNITS

! $Log$
! Revision 2.43  2014/04/22 00:05:34  vsnyder
! Remove unused USE name
!
! Revision 2.42  2014/03/20 01:39:08  vsnyder
! Add some tracing, use Value argument to search decl for base unit
!
! Revision 2.41  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.40  2013/10/09 23:39:23  vsnyder
! Put lit index in 'tree' field of unit declaration
!
! Revision 2.39  2013/09/19 23:25:28  vsnyder
! Add Base_Units
!
! Revision 2.38  2013/09/03 23:59:33  pwagner
! Added us (microsecond)
!
! Revision 2.37  2011/12/23 23:32:22  vsnyder
! Add MIFExtinctionv2 (as a comment)
!
! Revision 2.36  2011/11/11 00:30:59  vsnyder
! Add comments about L_ExtinctionV2 and L_MIFextinction
!
! Revision 2.35  2010/01/23 01:03:00  vsnyder
! Remove LogIWC
!
! Revision 2.34  2009/09/19 00:35:17  vsnyder
! Add phyq_logIceDensity and LotIWC
!
! Revision 2.33  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.32  2009/05/13 20:39:44  vsnyder
! Remove unused USE names
!
! Revision 2.31  2008/09/30 22:22:06  vsnyder
! Add iwc, gm3, mgm3
!
! Revision 2.30  2007/12/19 03:59:14  vsnyder
! Add a comment explaining why unused names are imported from Constants
!
! Revision 2.29  2007/12/06 20:37:57  vsnyder
! Move constants to Constants module to avoid sucking the parser infrastructure
! into off-line stuff that doesn't need it.
!
! Revision 2.28  2006/01/26 00:30:35  pwagner
! DU synonym for DobsonUnits
!
! Revision 2.27  2006/01/11 17:01:01  pwagner
! Made molcm2 default for colmabundance, DobsonUnits an alternate
!
! Revision 2.26  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.25  2004/06/03 22:58:20  vsnyder
! Add ms and kms units
!
! Revision 2.24  2004/03/26 01:31:24  vsnyder
! Add sqrtln2 constant
!
! Revision 2.23  2004/03/20 04:04:02  vsnyder
! Move Boltz and SpeedOfLight to physics
!
! Revision 2.22  2003/08/18 18:14:54  livesey
! Bug fix in declaration of iceDensity
!
! Revision 2.21  2003/08/16 00:34:02  vsnyder
! Use rad2deg instead of 180.0/Pi, push USE INTRINSIC down to procedure scope
!
! Revision 2.20  2003/01/26 04:41:25  livesey
! Added profiles
!
! Revision 2.19  2003/01/10 21:54:50  vsnyder
! Move SpeedOfLight from Geometry, put kinds on constants
!
! Revision 2.18  2003/01/07 23:43:44  livesey
! Added Gauss
!
! Revision 2.17  2002/12/02 23:00:19  vsnyder
! Add Sqrt Pi
!
! Revision 2.16  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.15  2002/09/27 00:27:52  vsnyder
! Remove Omega -- it's properly in Geometry.
! Move some USEs from module scope to procedure scope.
!
! Revision 2.14  2002/09/27 00:18:05  vsnyder
! Move USEs from module scope to procedure scope, add Boltzman constant
!
! Revision 2.13  2002/04/10 17:42:59  pwagner
! Added pctrhi unit
!
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
