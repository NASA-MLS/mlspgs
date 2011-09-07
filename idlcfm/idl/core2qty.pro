pro core2qty, core, qty, name, type, signal=signal, molecule=molecule, module=module, radiometer=radiometer
    case core.vertical of
        'P': begin
             vcoord = "zeta"
             surfs = double(-alog10(core.surfs))
             end
        'A': begin
             vcoord = "geodAltitude"
             surfs = double(1000 * core.surfs)
             end
        'N': begin
             vcoord = "none"
             surfs = double(core.surfs)
             end
        'G': begin
             vcoord = "gph"
             surfs = double(core.surfs)
             end
        'T': begin
             vcoor = "theta"
             surfs = double(core.surfs)
             end
        'X': begin
             vcoord = "integer"
             surfs = double(core.surfs)
             end
    endcase

    help, core

    tags = tag_names(core)
    nosurfs = long(core.nosurfs)
    noprofs = long(core.noprofs)
    noaux = long(core.noaux)
    if where(tags eq 'GEODANGLE') gt -1 then begin
        phi=double(temporary(core.geodangle))
        if core.coherence then phi = reform(phi, 1, noprofs) else phi = transpose(phi)
    endif
    if where(tags eq 'LOSANGLE') gt -1 then begin
        losangle=double(temporary(core.losangle))
        if core.coherence then losangle = reform(losangle, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = losangle
            losangle = s
        endelse
    endif
    if where(tags eq 'SZA') gt -1 then begin
        solarzenith=double(temporary(core.sza))
        if core.coherence then solarzenith = reform(solarzenith, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = solarzenith
            solarzenith = s
        endelse
    endif
    if where(tags eq 'LST') gt -1 then begin
        solarTime=double(temporary(core.lst))
        if core.coherence then solartime = reform(solartime, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = solartime
            solartime = s
        endelse
    endif
    if where(tags eq 'TIME') gt -1 then begin
        time=double(temporary(core.time))
        if core.coherence then time = reform(time, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = time
            time = s
        endelse
    endif
    if where(tags eq 'LAT') gt -1 then begin
        lat=double(temporary(core.lat))
        if core.coherence then lat = reform(lat, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = lat
            lat = s
        endelse
    endif
    if where(tags eq 'LON') gt -1 then begin
        lon=double(temporary(core.lon))
        if core.coherence then lon = reform(lon, 1, noprofs) else begin
            s = make_array(nosurfs, noprofs, /double)
            for i = 0, nosurfs-1 do s[i, *] = lon
            lon = s
        endelse
    endif

    if noAux gt 1 then begin
        module = core.auxinfo[0].module
        frequencies = core.auxinfo[*].v
        fcoord = "intermediatefrequency"
    endif
    ; The first reform is to turn value into a 2-dimensional array
    ; from maybe a 3-dimensional array or 1-dimensional array,
    ; so we can transpose.
    ; The sedond reform is to make sure value stays a 2-dimensional
    ; array after transpose.
    ; This is why I hate IDL. -haley
    value = reform(double(temporary(core.val)), noprofs, nosurfs * noaux)
    value = transpose(value)
    value = reform(value, nosurfs * noaux, noprofs)

    if n_elements(surfs) ge 0 then begin
        if core.coherence then begin
            surfs = reform(surfs, nosurfs, 1)
        endif else surfs = transpose(surfs) 
        ; because input data array shape is inconsistent
        ; at this point, if surfs has a trailing dimension of 1
        ; it may or may not be dropped depending on which if
        ; clause it passes through above. However, as long as
        ; the number of elements is enough, the receiving
        ; side will force the elements into the correct shape
        ; on the server side.
    endif

    createquantity, qty, name, type, 0L, value=value, logbasis=logbasis, $
    startval=startval, noChans=noAux, surfs=surfs, verticalCoordinate=vcoord, $
    phi=phi, losangle=losangle, longitude=lon, geodlat=lat, time=time, $
    noProfs=noProfs, noSurfs=noSurfs, signal=signal, $
    solarzenith=solarzenith, solartime=solartime, module=module, $
    frequencies=frequencies, frequencyCoordinate=fcoord, $
    badval=double(core.baddata), molecule=molecule, radiometer=radiometer
end
