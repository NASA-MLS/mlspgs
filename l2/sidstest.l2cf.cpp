;@G
; $Id$

; This is a level 2 configuration file for an interim v0.5 version of
; the EOS MLS Level 2 software.

; This intent of this l2cf is simply to provide a sids style radiance
; dataset for band 2.

; -----------------------------------------------------  MLSSignals  -----
BEGIN MLSSignals

#include "mls_signals.l2cf"

END MLSSignals

; ---------------------------------------------------  Spectroscopy  -----
BEGIN Spectroscopy

#include "big_spect_cat.l2cf"

END Spectroscopy

; -------------------------------------------------  GlobalSettings  -----
BEGIN GlobalSettings

    ;; In reality we probably don't need any of these.
    VersionComment="L2CF Version 0-1-0 for version 0.1 of the MLS software"
    InputVersionString="v0-1-0-test"
    OutputVersionString="v0-1-0-test"
    AllowClimatologyOverloads=FALSE

    ;; Construct vGrids

    fwmIntegrationGrid: vGrid, coordinate=Zeta, type=Logarithmic, $
      start = 1000mb, formula=[ 65:48, 40:24, 20:12]
    ;; Forward model tangent grid needs points below the surface.
    fwmTangentGrid: vGrid, coordinate=Zeta, type=Logarithmic, $
      start = 10.0e9 mb, formula=[8:1, 32:24, 20:12, 10:6]

    vGridStandard: vGrid, coordinate=Zeta, type=Logarithmic, $
                   start=1000mb, formula=[25:12, 24:6]
    vGridRefGPH: vGrid, coordinate=Zeta, type=Explicit, values=100mb

;; Stuff to start up the forward model from Zvi's and Bill's files
    l2load,   zvi='/user5/zvi/n_seez'
;; End of stuff to start up the forward model from Zvi's and Bill's
;; files

    forwardModelGlobal, $
      antennaPatterns='/user5/vsnyder/mlspgs/l2/aaap.umls', $
      filterShapes='/user5/vsnyder/mlspgs/l2/normalized_filter_banks.dat', $
      pointingGrids='/user1/livesey/mlspgs/l2/full_r2_b234_lsb'

    testFwm: forwardModel, type=full, signals = [sig_band2], $
      molecules = [h2o], frequency=182610.091 MHz, channels = [1], $
      integrationGrid=fwmIntegrationGrid, tangentGrid=fwmTangentGrid, $
      /do_conv, /do_freq_avg

END GlobalSettings


; ----------------------------------------------------  ReadApriori  -----
BEGIN ReadApriori

    ;; Currently this reads apriori data from l2gp files, designed for
    ;; 1 January 2000.  Later we should use gridded data, or at least
    ;; stuff for 1996

    l2gpTempInput: l2gp, file='/nas/bigdata/livesey/misc/temp_input.l2gp', $
      swath='Temperature'
    l2gpH2OInput: l2gp, file='/nas/bigdata/livesey/misc/h2o_input.l2gp', $
      swath='H2O'
    l2gpREFGPHInput: l2gp, file='/nas/bigdata/livesey/misc/refgph_input.l2gp', $
      swath='refGPH'

;    gridTemp: gridded, file='/nas/bigdata/livesey/misc/temp_grid.dat', $
;      field='Temperature', origin=climatology
    
END ReadApriori


; ---------------------------------------------------  MergeApriori  -----
BEGIN MergeApriori
    ;; Nothing to do here
END MergeApriori


; ----------------------------------------------------  ChunkDivide  -----
BEGIN ChunkDivide
    ;; Just a pretty typical setup, except no overlaps to avoid bug

 ;; IdealLength =             0.00416666667  orbits
  ;; IdealLength =             0.25 orbits
    IdealLength =             0.03  orbits ; Chosen to give 5MAFs in first chunk
    Overlap =                 0 MAFs

    HomeGeodAngle =           1.5 degrees ; Chosen to give 5MAFs in the first chunk
    HomeModule =              GHz

    ScanLowerLimit =          -10 km : 20 km
    ScanUpperLimit =          40 km : 100 km

    CriticalScanningModules = Both
    CriticalBands =           "R1A:118.B1:PT or R1B:118.B21:PT and R2:190.B*"

    MaxGap =                  3 minutes
END ChunkDivide


; ------------------------------------------------------  Construct  -----
BEGIN Construct

    ;; First construct hGrids
    hGridStandard: hGrid, module=THz, type=height, height=15km

    ;; Now some vector quantities, main state vector stuff first.
    temp: quantity, vGrid=vGridStandard, hGrid=hGridStandard, type=temperature
    refGPH: quantity, vGrid=vGridRefGPH, hGrid=hGridStandard, type=refGPH
    H2O: quantity, vGrid=vGridStandard, hGrid=hGridStandard, $
       type=vmr, molecule=H2O, unit=ppmv
    ptanGHz: Quantity, type=ptan, module=GHz

    ;; Now the more boring state stuff
    b2elevOffset: Quantity, type = elevOffset, signal = sig_band2
    earthReflectivity: quantity, type = earthRefl
    orbitIncline: Quantity, type = orbitIncline
    scGeocAlt: Quantity, type = scGeocAlt, module= SC
    spaceRadiance: quantity, type = spaceRadiance
    scECI: quantity, type=scECI, module=SC
    scVel: quantity, type=scVel, module=SC
    tpECIGHz: quantity, type=tngtECI, module=GHz
    losVelGHz: quantity, type=losVel, module=GHz

    ;; Now some measurement quantities
    band2: Quantity, signal=sig_band2, type=radiance
    tpGeocAltGHz: Quantity, type=tngtGeocAlt, module=GHz

    ;; Now define a state vector in two parts main an extra
    stateTemplate: VectorTemplate, quantities = $
      [temp, refGPH, ptanGHz,H2O]
    extraStateTemplate: VectorTemplate, quantities =  $
      [b2elevOffset, earthReflectivity, orbitIncline, scGeocAlt,  $
       spaceRadiance, scECI, scVel, losVelGHz, tpECIGHz]

    ;; Now define templates for the measurement vectors.
    radianceTemplate: VectorTemplate, quantities = [band2]
    miscMeasTemplate: VectorTemplate, quantities = [tpGeocAltGHz]
    

END Construct


; ------------------------------------------------------------  Fill -----
BEGIN Fill

    ;; First create our vectors
    state:      Vector, template = stateTemplate
    extraState: Vector, template = extraStateTemplate

    radiance: Vector, template = radianceTemplate
    miscMeas: Vector, template = miscMeasTemplate

    ;; Fill most of the state stuff from l2gp
    Fill, quantity = state.temp, method=l2gp, sourcel2gp=l2gpTempInput
    Fill, quantity = state.refGPH, method=l2gp, sourceL2GP=l2gpRefGPHInput
    Fill, quantity = state.h2o, method=l2gp, sourceL2GP=l2gpH2OInput

    ;; Fill the miscMeas stuff from l1b, we need this to guess ptan
    Fill, quantity = miscMeas.tpGeocAltGHz, method = l1b

    ;; Now use a hydrostatic calculation to fill ptanGHz
    Fill, quantity=state.ptanGHz, method=hydrostatic, $
      h2oQuantity=state.h2o, $
      temperatureQuantity=state.temp, $
      refGPHQuantity=state.refGPH, $
      geocAltitudeQuantity=miscMeas.tpGeocAltGHz, $
      maxIterations=4

    ;; Now fill the extra state stuff, mostly explicitly
    Fill, quantity = extraState.b2elevOffset, $
      method = explicit, explicitValues = 0 degrees
;    Snoop, comment="Towards the start of fill"

    Fill, quantity = extraState.earthReflectivity, $
      method = explicit, explicitValues = 0.05 ; CHECK THIS NUMBER !!!!!
    Fill, quantity = extraState.orbitIncline, $
      method = explicit, explicitValues = 98.145 degrees
;    Snoop, comment="In the middle of fill"

    Fill, quantity = extraState.scGeocAlt, method = l1b
    Fill, quantity = extraState.scECI, method = l1b
    Fill, quantity = extraState.scVel, method = l1b
    Fill, quantity = extraState.tpECIGHz, method = l1b

    Fill, quantity = extraState.losVelGHz, tngtECI=extraState.tpECIGHz, $
      scECI=extraState.scECI, scVel=extraState.scVel, method=special

    Fill, quantity = extraState.spaceRadiance, $
      method =explicit, explicitValues = 2.735 K

; Snoop, comment="At the end of fill"
    
END Fill


; -------------------------------------------------------  Retrieve  -----

BEGIN Retrieve

    sids, fwdModelIn = state, fwdModelExtra = extraState, fwdModelOut = radiance, $
      forwardModel=testFwm

END Retrieve


; -----------------------------------------------------------  Join  -----
BEGIN Join

    tempL2gp: l2gp, source=state.temp, swath='Temperature'
    h2oL2gp: l2gp, source=state.h2o, swath='H2O'
    refgphL2GP: l2gp, source=state.refgph, swath='refGPH'

    ptanGHzl2aux: l2aux, source=state.ptanGHz, sdName='ptanGHz'
    tpECIl2aux: l2aux, source=extrastate.tpECIGHz, sdName='tpECIGHz'
    losVelL2aux: l2aux, source=extrastate.losVelGHz, sdName='losVelGHz'
    b2l2aux: l2aux, source = radiance.band2, sdName='R2:240.B2F:H2O.S0.FB25-2'

END Join


; ---------------------------------------------------------  Output  -----
BEGIN Output

    output, file='l2gp_temp', quantities=tempL2gp,type=l2gp
    output, file='l2gp_gph', quantities=refgphL2gp,type=l2gp
    output, file='l2gp_h2o', quantities=h2oL2gp,type=l2gp
    
    output, file='l2aux_full', quantities= [ ptanGHzl2aux, b2l2aux, $
                                            tpECIl2aux, losVell2Aux], type=l2aux

END Output

; ==========================================================================

; $Log$
; Revision 2.1  2001/04/05 02:08:17  vsnyder
; Initial commit
;
