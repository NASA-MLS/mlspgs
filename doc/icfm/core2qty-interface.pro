pro core2qty, core, qty, name, type, signal=signal, molecule=molecule, module=module
    ; core: AtGod quantity, usually created by AtGod_CreateCore function
    ; qty: output, ICFM Quantity data type that is equivalent to core
    ; name: unlike core, qty needs a name
    ; type: unlike core, qty also needs a type. Please refer to Appendix 
    ; for applicable types.
    ; signal (optional): some quantities have an associated signal. Please
    ; refer to Appendix for applicable signals.
    ; molecule (optional): some quantities have an associated molecule.
    ; Please refer to Appendix for applicable molecules.
    ; module (optional): some quantities are associated with a module. 
    ; Module is one of "GHz", "THz", and "sc".
end
