pro CreateQuantity, quantity, name, type, instanceOffset, value=value, $
logbasis=logbasis, startValue=startValue, badValue=badValue, $
molecule=molecule, module=instrumentModule, signal=signal, $
radiometer=radiometer, frequencies=frequencies, $
frequencyCoordinate=frequencyCoordinate, noChans=noChans, surfs=surfaces, $
nosurfs=noSurfs, verticalCoordinate=verticalCoordinate, noprofs=noprofs, $
phi=phi, geodlat=geodlat, longitude=lon, losAngle=losAngle, $
solarZenith=solarZenith, solarTime=solarTime, time=time, mask=mask, $
coherence=coherence
    ; quantity: output, the quantity to be returned
    ; name: name of the quantity
    ; type: a string. Please refer to Appendix for applicable types.
    ; instanceOffset: index of the first non-overlapped instance
    ; value (optional): the value for this quantity
    ; logbasis (optional): whether the value is on a log scale. Default is false.
    ; startValue (optional): only valid when logbasis is true. 
    ; badValue (optional): the constant representing an invalid value point 
    ; in the value array
    ; molecule (optional): the molecule associated with this quantity if any.
    ; Please refer to Appendix for a list of applicable molecules.
    ; module (optional): the instrument module associated to this quantity
    ; if any. Module is a string, one of "GHz", "THz", "sc".
    ; signal (optional): the signal associated with this quantity if any.
    ; Please refer to Appendix for a list of applicable signals.
    ; radiometer (optional): if signal is not needed, there can still a 
    ; a radiometer associated with this quantity. Please refer to the 
    ; signal list in the Appendix for applicable radiometers.
    ; frequencies (optional): 
    ; frequencyCoordinate (optional): valid only when frequencies are present
    ; noChans (optional): needed only when frequencies are not present. When
    ; frequencies are present, noChans is set to the size of frequencies. 
    ; Otherwise, this value is 1.
    ; surfs (optional): associated pressure surfaces
    ; nosurfs (optional): number of pressure surfaces
    ; verticalCoordinate (optional): only valid when surfs is present.
    ; Please refer to Appendix for a list of valid vertical coordinate.
    ; noProfs (optional): a.k.a. noInstances
    ; phi (optional)
    ; geodlat (optional):
    ; lon (optional):
    ; losAngle (optional):
    ; solarZenith (optional):
    ; solarTime (optional): 
    ; time (optional):
    ; coherence (optional): default is true
end 
