; Just a mockup of how one might use the idl-callable forward model
; ---------------------------------- -------------------------------------------
pro mockup
l2aux = '/data/emls/l2aux/v03.30/2009/040/MLS-Aura_L2AUX-DGM_v03-30-c01_2009d040.h5'
l1boa = '/data/emls/l1b/v03.30/2009/040/MLS-Aura_L1BOA_v03-30-c01_2009d040.h5'
l1brad = '/data/emls/l1b/v03.30/2009/040/MLS-Aura_L1BRADG_v03-30-c01_2009d040.h5'
filenames = [l1boa, l1brad, l2aux]
maf = 214
band = 9
ptan = reademlsmifq(filenames=filenames, 'GHz.ptan-Core', vertical='A', nomafs=1, $
firstmaf=maf)
createvector, stateExtra, 'stateExtra'
core2qty, ptan, temp, 'ptan', "ptan", module="GHz"
addquantity2vector, stateExtra, temp
geodalt = reademlsmifq(filenames=filenames, 'GHz/GeodAlt', firstmaf=maf, nomafs=1, $
vertical='N')
core2qty, geodalt, temp, 'geodaltitude', 'tngtgeodalt', module='GHz'
addquantity2vector, stateExtra, temp
losvelghz = reademlsmifq(filenames=filenames, 'GHz/LosVel', firstmaf=maf, nomafs=1, $
vertical='N')
core2qty, losvelghz, temp, 'losVelGHz', 'LOSVel', module='GHz'
addquantity2vector, stateExtra, temp
orbincl = reademlsmifq(filenames=filenames, 'sc/OrbIncl', firstmaf=maf, nomafs=1, $
vertical='N', module='sc', /limitedgeolocation)
core2qty, orbincl, temp, 'orbincl', 'orbitinclination', module='sc'
addquantity2vector, stateextra, temp
phitan = ptan
phitan.val = ptan.geodangle
core2qty, phitan, temp, 'phitan', 'phitan', module='GHz'
addquantity2vector, stateextra, temp
spaceRadiance = AtGod_CreateCore ( vertical='N', noProfs=1 )
spaceRadiance.val = 2.735
core2qty, spaceradiance, temp, 'spacerad', 'spaceradiance'
addquantity2vector, stateextra, temp
earthReflectivity = AtGod_CreateCore ( vertical='N', noProfs=1 )
earthReflectivity.val = 0.05
core2qty, earthreflectivity, temp, 'earthrefl', 'earthRefl'
addquantity2vector, stateextra, temp
scgeocalt = reademlsmifq (filenames=filenames, 'sc/GeocAlt', vertical='N', $
firstmaf=maf, nomafs=1, module='sc', /limitedgeolocation)
core2qty, scgeocalt, temp, 'scGeocAlt', 'scgeocalt', module='sc'
addquantity2vector, stateextra, temp
ghzgeocalt = reademlsmifq (filenames=filenames, 'GHz/GeocAlt', vertical='N', $
firstmaf=maf, nomafs=1, module='GHz')
core2qty, ghzgeocalt, temp, 'geocAlt', 'tngtgeocalt', module='GHz'
addquantity2vector, stateextra, temp
sbfFile = '/data/emls/l2cal/MLS-Aura_L2Cal-SBFraction_v3-0-2_0000d000.l2cf'
band='B9F'
band9 = reademlsmifq(filenames=filenames, signal=band, vertical='A', nomafs=1, $
firstmaf=maf)
createvector, measurement, 'measurement'
core2qty, band9, temp, 'band9', 'radiance', signal='R3:240.B9F:CO'
addquantity2vector, measurement, temp
sbf9l = atgod_createcore(vertical='N', noProfs=1, auxinfo=band9.auxinfo)
sbf9l.val = readl2cfarray(sbffile, $
prefix='!define(limbradsbfraction9L,{', suffix='})!dnl', extrasuffix=',/spread')
core2qty, sbf9l, temp, 'lsf9l', 'limbsidebandfraction', $
signal='R3:240.B9LF:CO.S0.FB25-9'
addquantity2vector, stateextra, temp
sbf9u = atgod_createcore(vertical='N', noProfs=1, auxinfo=band9.auxinfo)
sbf9u.val = readl2cfarray(sbffile, $
prefix='!define(limbradsbfraction9U,{', suffix='})!dnl', extrasuffix=',/spread')
core2qty, sbf9u, temp, 'lsf9u', 'limbsidebandfraction', signal='R3:240.B9UF:CO'
addquantity2vector, stateextra, temp
elevoffset9l = atgod_createcore(vertical='N', noProfs=1, auxinfo=band9.auxinfo)
elevoffset9l.val = 0
core2qty, elevoffset9l, temp, 'eo9l', 'elevoffset', signal='R3:240.B9LF:CO'
addquantity2vector, stateextra, temp
elevoffset9u = atgod_createcore(vertical='N', noprofs=1, auxinfo=band9.auxinfo)
elevoffset9u.val = 0
core2qty, elevoffset9u, temp, 'eo9u', 'elevoffset', signal='R3:240.B9UF:CO'
addquantity2vector, stateextra, temp
date = '2009d040'
version = 'v03.30'
CO = ReadEMLSDaySeries ( date, date, version=version, order='cycle desc', $
    product='CO', /quiet )
;; Select only the first 'noProfs' profiles
noProfs = maf + 6
CO = ExtractProfiles ( CO,  noProfs)
core2qty, co, temp, 'CO', 'vmr', molecule='CO'
createvector, state, 'state'
addquantity2vector, state, temp
temperature = reademlsdayseries(date, date, version=version, order='cycle desc', $
product='Temperature', /quiet)
temperature = extractprofiles(temperature, noprofs)
core2qty, temperature, temp, 'Temperature', 'temperature'
addquantity2vector, stateextra, temp
gph = reademlsdayseries(date, date, version=version, order='cycle desc', $
product='GPH', /quiet)
gph = extractprofiles(gph, noprofs)
core2qty, gph, temp, 'GPH', 'gph'
addquantity2vector, stateextra, temp
refgph = atgod_createcore(vertical='P', nosurfs=1, noprofs=temperature.noprofs)
refgph.surfs = 100.0D
refgph.val = 16000.0D
refgph.lat = temperature.lat
refgph.lon = temperature.lon
refgph.sza = temperature.sza
refgph.day = temperature.day
refgph.time = temperature.time
refgph.lst = temperature.lst
refgph = create_struct(['losangle', 'geodangle'], temperature.losangle, $
temperature.geodangle, refgph)
core2qty, refgph, temp, 'refGPH', 'refGPH'
addquantity2vector, stateextra, temp
common pvmlib, config
if n_elements(config) eq 0 then Initpvmlib
msgtag=200L
signalFile='/users/honghanh/mlspgs/cfm/signals.l2cf'
configFile='/users/honghanh/mlspgs/idlcfm/config.l2cf'
mafNo = 0L
spectroscopy='/data/emls/l2cal/MLS-Aura_L2Cal-Spectroscopy-PFA_v3-0-4_0000d000.h5'
antennaPatterns='/data/emls/l2cal/MLS-Aura_L2Cal-AAAP_v2-0-0_0000d000.txt'
filterShapes='/data/emls/l2cal/MLS-Aura_L2Cal-Filters_v3-0-2_0000d000.txt'
DACSFilterShapes='/data/emls/l2cal/MLS-Aura_L2Cal-DACSFilters_v1-5-1_0000d000.txt'
pointingGrids='/data/emls/l2cal/MLS-Aura_L2Cal-PFG_v3-0-0_0000d000.txt'
pfa=['/data/emls/l2cal/PFA_R5V_FS-04.h5', $
'/data/emls/l2cal/PFA_R5H_FS-04.h5', $
'/data/emls/l2cal/PFA_R4_FS-04.h5', $
'/data/emls/l2cal/PFA_R3_FS-04.h5', $
'/data/emls/l2cal/PFA_R2_FS-04.h5', $
'/data/emls/l2cal/PFA_R1B_FS-04.h5', $
'/data/emls/l2cal/PFA_R1A_FS-04.h5', $
'/data/emls/l2cal/PFA_DACS_FS-04.h5']
hdf5l2pc= $
['/data/emls/l2cal/l2pc_30H6/MLS-Aura_L2Cal-L2PC-band9-LATSCALARHIRESO3HR_v3-00-HO-06_m02.h5']

; ------------------------------ tid -------------------------------------------
; tid is the task id the server prints when it starts up
; it also saves it to a text file named tid.txt
; uncomment and use one of the following lines to set tid for this run
; tid = 262146 ; 1 ; replace this with the tid that server prints when it starts
; or else use geetenv, like the following
; tid = getenv('tid')
; ---------------------------------- -------------------------------------------
icfm_setup, tid,info, signalfile, configfile, spectroscopy, $
antennaPatterns, filterShapes, dacsFilterShapes, pointingGrids, $
pfa, hdf5l2pc
forwardmodel, tid, info, mafNo, state, stateExtra, measurement, $
output
icfm_cleanup, tid, info
end
; ---------------------------------- -------------------------------------------
; $Log$
