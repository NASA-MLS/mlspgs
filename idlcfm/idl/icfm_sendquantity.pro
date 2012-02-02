pro ICFM_SendQuantity, qty, tid=tid, info=info, msgtag=msgtag, packonly=packonly
    if size(qty, /type) ne 8 then MyMessage, /error, "qty must be created via createquantity procedure"

    tags = tag_names(qty)
    has_name = where(tags eq 'NAME') gt -1
    has_type = where(tags eq 'TYPE') gt -1
    has_offset = where(tags eq 'INSTANCEOFFSET') gt -1
    has_logbasis = where(tags eq 'LOGBASIS') gt -1
    has_startval = where(tags eq 'STARTVALUE') gt -1
    has_badval = where(tags eq 'BADVALUE') gt -1
    has_molecule = where(tags eq 'MOLECULE') gt -1
    has_module = where(tags eq 'MODULE') gt -1
    has_signal = where(tags eq 'SIGNAL') gt -1
    has_radiometer = where(tags eq 'RADIOMETER') gt -1
    has_frequencies = where(tags eq 'FREQUENCIES') gt -1
    has_fCoord = where(tags eq 'FREQUENCYCOORDINATE') gt -1
    has_noChans = where (tags eq 'NOCHANS') gt -1
    has_surfs = where (tags eq 'SURFACES') gt -1
    has_vCoor = where (tags eq 'VERTICALCOORDINATE') gt -1
    has_noSurfs = where(tags eq 'NOSURFS') gt -1
    has_phi = where (tags eq 'PHI') gt -1
    has_geodlat = where(tags eq 'GEODLAT') gt -1
    has_lon = where (tags eq 'LONGITUDE') gt -1
    has_losAngle = where (tags eq 'LOSANGLE') gt -1
    has_solarZenith = where(tags eq 'SOLARZENITH') gt -1
    has_solarTime = where(tags eq 'SOLARTIME') gt -1
    has_time = where(tags eq 'TIME') gt -1
    has_noprofs = where(tags eq 'NOINSTANCES') gt -1
    has_instancelen = where(tags eq 'INSTANCELEN') gt -1
    has_value = where(tags eq 'VALUE') gt -1
    has_mask = where(tags eq 'MASK') gt -1
    has_coherent = where (tags eq 'COHERENT') gt -1
    has_stacked = where (tags eq 'STACKED') gt -1

    l29 = [has_name, has_type, has_offset, has_coherent, $
           has_stacked, $
           has_logbasis, has_startval, has_badval, $
           has_molecule, has_module, has_signal, $
           has_radiometer, has_frequencies, has_fCoord, $
           has_noChans, has_surfs, has_vCoor, has_noSurfs, $
           has_phi, has_geodlat, has_lon, has_losAngle, $
           has_solarzenith, has_solartime, has_time, $
           has_noprofs, has_instancelen, has_value, has_mask]
    l29 = long(temporary(l29)) and 'ffff'xL
    pvm_pack_idltypeinfo, l29
    pvm_pack_idlvariable, l29

    if has_name then begin
        pvm_pack_idltypeinfo, qty.name
        pvm_pack_idlvariable, qty.name
    endif

    if has_type then begin
        pvm_pack_idltypeinfo, qty.type
        pvm_pack_idlvariable, qty.type
    endif else MyMessage, /error, "qty doesn't have 'type' attribute"

    if has_offset then begin
        pvm_pack_idltypeinfo, qty.instanceOffset
        pvm_pack_idlvariable, qty.instanceOffset
    endif else MyMessage, /error, "qty doesn't have 'instanceOffset' attribute"

    if has_coherent then begin
        pvm_pack_idltypeinfo, qty.coherent
        pvm_pack_idlvariable, qty.coherent
    endif else MyMessage, /error, "qty doesn't have 'coherent' attribute"

    if has_stacked then begin
        pvm_pack_idltypeinfo, qty.stacked
        pvm_pack_idlvariable, qty.stacked
    endif else MyMesssage, /error, "qty doesn't have 'stack' attribute"

    if has_logbasis then begin
        pvm_pack_idltypeinfo, qty.logbasis
        pvm_pack_idlvariable, qty.logbasis
    endif

    if has_startval then begin
        pvm_pack_idltypeinfo, qty.startValue
        pvm_pack_idlvariable, qty.startValue
    endif

    if has_badval then begin
        pvm_pack_idltypeinfo, qty.badValue
        pvm_pack_idlvariable, qty.badValue
    endif

    if has_molecule then begin
        pvm_pack_idltypeinfo, qty.molecule
        pvm_pack_idlvariable, qty.molecule
    endif

    if has_module then begin
        pvm_pack_idltypeinfo, qty.module
        pvm_pack_idlvariable, qty.module
    endif

    if has_signal then begin
        pvm_pack_idltypeinfo, qty.signal
        pvm_pack_idlvariable, qty.signal
    endif

    if has_radiometer then begin
        pvm_pack_idltypeinfo, qty.radiometer
        pvm_pack_idlvariable, qty.radiometer
    endif

    if has_nochans then begin
        pvm_pack_idltypeinfo, qty.noChans
        pvm_pack_idlvariable, qty.noChans
    endif

    if has_fcoord then begin
        pvm_pack_idltypeinfo, qty.frequencyCoordinate
        pvm_pack_idlvariable, qty.frequencyCoordinate
    endif

    if has_frequencies and has_nochans and has_fcoord then begin
        pvm_pack_idltypeinfo, qty.frequencies
        pvm_pack_idlvariable, qty.frequencies
    endif else if has_frequencies then MyMessage, /error, "frequencies must go with nochans and frequencyCoordinate"

    if has_nosurfs then begin
        pvm_pack_idltypeinfo, qty.nosurfs
        pvm_pack_idlvariable, qty.nosurfs
    endif

    if has_noprofs then begin
        pvm_pack_idltypeinfo, qty.noinstances
        pvm_pack_idlvariable, qty.noinstances
    endif

    if has_phi and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.phi
        pvm_pack_idlvariable, qty.phi
    endif else if has_phi then MyMessage, /error, "phi must be accompanied by noInstances"

    if has_geodlat and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.geodlat
        pvm_pack_idlvariable, qty.geodlat
    endif else if has_geodlat then MyMessage, /error, "geodlat must be accompanied by noInstances"

    if has_lon and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.longitude
        pvm_pack_idlvariable, qty.longitude
    endif else if has_longitude then MyMessage, /error, "longitude must be accompanied by noInstances"

    if has_losangle and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.losangle
        pvm_pack_idlvariable, qty.losangle
    endif else if has_losangle then MyMessage, /error, "losangle must be accompanied by noInstances"

    if has_solarzenith and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.solarzenith
        pvm_pack_idlvariable, qty.solarzenith
    endif else if has_solarzenith then MyMessage, /error, "solarzenith must be accompanied by noInstances"

    if has_solartime and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.solartime
        pvm_pack_idlvariable, qty.solartime
    endif else if has_solartime then MyMessage, /error, "solartime must be accompanied by noInstances"

    if has_time and has_noprofs then begin
        pvm_pack_idltypeinfo, qty.time
        pvm_pack_idlvariable, qty.time
    endif else if has_time then MyMessage, /error, "time must be accompanied by noInstances"

    if has_vcoor then begin
        pvm_pack_idltypeinfo, qty.verticalCoordinate
        pvm_pack_idlvariable, qty.verticalCoordinate
    endif

    if has_surfs and has_nosurfs and has_vcoor then begin
        pvm_pack_idltypeinfo, qty.surfaces
        pvm_pack_idlvariable, qty.surfaces
    endif else if has_surfs then MyMessage, /error, "surfs must be accompanied by nosurfs and verticalCoordinate"

    if has_instancelen then begin
        pvm_pack_idltypeinfo, qty.instanceLen
        pvm_pack_idlvariable, qty.instanceLen
    endif

    if has_value then begin
        pvm_pack_idltypeinfo, qty.value
        pvm_pack_idlvariable, qty.value
    endif

    if has_mask and has_value then begin
        pvm_pack_idltypeinfo, qty.mask
        pvm_pack_idlvariable, qty.mask
    endif else if has_mask then MyMessage, /error, "mask must be accompanied by value"

    if n_elements(packonly) eq 0 then begin
        ; send data
        if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
        else info = pvm_mcast(tids=tid, msgtag=msgtag)
    endif


end
