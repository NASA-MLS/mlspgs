; This subroutine utilize some subroutines in AtGod
pro ICFM_Setup, tid, info, signalfile, configfile, spectroscopy, antennaPatterns, filterShapes, dacsFilterShapes, pointingGrids, pfa, hdf5l2pc
    
    if n_elements(signalfile) eq 0 then MyMessage, /error, "Missing signalfile"
    if n_elements(configfile) eq 0 then MyMessage, /error, "Missing configfile"
    if n_elements(spectroscopy) eq 0 then MyMessage, /error, "Missing spectroscopy"
    if n_elements(antennaPatterns) eq 0 then MyMessage, /error, "Missing antennaPatterns"
    if n_elements(filterShapes) eq 0 then MyMessage, /error, "Missing filterShapes"
    if n_elements(dacsFilterShapes) eq 0 then MyMessage, /error, "Missing dacsFilterShapes"
    if n_elements(pointingGrids) eq 0 then MyMessage, /error, "Missing pointingGrids"
    if n_elements(pfa) eq 0 then MyMessage, /error, "Missing pfa"
    if n_elements(hdf5l2pc) eq 0 then MyMessage, /error, "Missing hdf5l2pc"

    bufid = pvm_initsend()

    pvm_pack_idltypeinfo, 0L ;sig_setup
    pvm_pack_idlvariable, 0L

    pvm_pack_idltypeinfo, signalfile
    pvm_pack_idlvariable, signalfile

    pvm_pack_idltypeinfo, configfile
    pvm_pack_idlvariable, configfile

    pvm_pack_idltypeinfo, spectroscopy
    pvm_pack_idlvariable, spectroscopy

    pvm_pack_idltypeinfo, antennaPatterns
    pvm_pack_idlvariable, antennaPatterns

    pvm_pack_idltypeinfo, filterShapes
    pvm_pack_idlvariable, filterShapes

    pvm_pack_idltypeinfo, dacsFilterShapes
    pvm_pack_idlvariable, dacsFilterShapes

    pvm_pack_idltypeinfo, pointingGrids
    pvm_pack_idlvariable, pointingGrids

    a = long(n_elements(pfa))
    pvm_pack_idltypeinfo, a
    pvm_pack_idlvariable, a
    for i=0,a-1 do begin
        pvm_pack_idltypeinfo, pfa[i]
        pvm_pack_idlvariable, pfa[i]
    endfor

    a = long(n_elements(hdf5l2pc))
    pvm_pack_idltypeinfo, a
    pvm_pack_idlvariable, a
    for i=0,a-1 do begin
        pvm_pack_idltypeinfo, hdf5l2pc[i]
        pvm_pack_idlvariable, hdf5l2pc[i]
    endfor

    msgtag = 200L
    if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
    else info = pvm_mcast(tids=tid, msgtag=msgtag)
end
