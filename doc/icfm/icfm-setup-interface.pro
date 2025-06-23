pro ICFM_Setup, tid, info, signalfile, configfile, spectroscopy, $
antennaPatterns, filterShapes, dacsFilterShapes, pointingGrids, $
pfa, hdf5l2pc
    ; This procedures allocate common databases of different forward
    ; model calls on the server.
    ; tid: server's tid
    ; info: output, status report info
    ; signalfile: absolute file name of the signal file
    ; configfile: absolute file name
    ; spectroscopy: absolute file name of spectroscopy file
    ; antennapatterns: absolute file name of the antenna file
    ; filterShapes: absolute file name
    ; dacsFilterShapes: absolute file name
    ; pointingGrids: absolute file name
    ; pfa: an array of absolute file names
    ; hdf5l2pc: an array of absolute file names
end
