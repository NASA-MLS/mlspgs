pro ForwardModel, tid, info, mafNo, state, stateExtra, measurement, output
    ; tid: server's tid
    ; info: output, status message of the call to foward model
    ; mafNo: the maf number of the data to compute
    ; state: input vector
    ; stateExtra: input vector
    ; measurement: input vector, but no values needed for quantities of this
    ; output: variable to hold output data. Output will have as many 
    ; quantities as measurement. However, those quantities are AtGod core
    ; data type.
end
